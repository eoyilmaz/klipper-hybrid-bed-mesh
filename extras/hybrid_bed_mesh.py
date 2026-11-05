"""Hybrid Eddy rapid scan + tap bed mesh.

An Eddy probe mounted behind the nozzle (positive ``y_offset``) can not reach
the front strip of the bed, the coil runs out of Y travel before the nozzle
does. ``BED_MESH_HYBRID_CALIBRATE`` builds a single bed mesh out of two
passes:

1. The eddy-reachable part of the mesh is rapid scanned with
   ``BED_MESH_CALIBRATE METHOD=rapid_scan``.
2. The front rows, plus the first scanned row as an overlap row, are probed
   with ``METHOD=tap``, which uses the nozzle itself and so has no XY offset.

The tap rows are corrected by the difference between the tap and the scan
values on the overlap row, so the front strip lines up with the scanned part,
and then both are merged into one evenly spaced mesh. The correction is either
the average difference (``offset``) or a straight line fitted along X
(``linear``), the latter also takes out an X slope between the nozzle and the
eddy coil that axis twist compensation leaves behind, like a gantry twist that
changes along Y.

The pure planning / merging functions have no Klipper dependencies, so they
can be unit tested standalone.
"""

# Standard Library Imports
from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

# Klipper Imports
from extras.bed_mesh import BedMeshError, ZMesh
from extras.probe import ProbePointsHelper

if TYPE_CHECKING:
    from collections.abc import Sequence

    from configfile import ConfigWrapper
    from gcode import GCodeCommand
    from klippy import Printer

    from extras.bed_mesh import BedMesh

logger = logging.getLogger(__name__)

# BED_MESH_HYBRID_CALIBRATE parameters handed over to the rapid scan pass
SCAN_PASS_THROUGH = ("ALGORITHM", "MESH_PPS", "HORIZONTAL_MOVE_Z", "SCAN_SPEED")

# the start of the eddy probe's error when it can't find the tap in the data,
# only these taps are retried
TAP_ERROR = "Unable to detect tap"

CALIBRATE_HELP = (
    "Rapid scan the bed mesh with the eddy probe and tap the front strip it "
    "can not reach"
)

# bed_mesh needs at least this many rows for bicubic interpolation
MIN_SCAN_ROWS = 4

# how the tap rows are lined up with the scan on the overlap row
OVERLAP_ALIGNMENTS = ("linear", "offset")


@dataclass
class HybridPlan:
    """The Y rows of a hybrid mesh and how they are split between the passes.

    Attributes:
        rows (list[float]): The Y coordinates of all mesh rows, front to back.
        split (int): The index of the overlap row, the first row the eddy
            coil can reach. Rows ``0..split`` are tapped, rows
            ``split..len(rows)-1`` are scanned.
        step (int): The scan pass samples ``step`` rows per mesh row, so that
            it has at least :data:`MIN_SCAN_ROWS` rows, only every ``step``-th
            scanned row ends up in the mesh.
    """

    rows: list[float]
    split: int
    step: int

    @property
    def scan_row_count(self) -> int:
        """int: The number of rows the rapid scan pass probes."""
        return (len(self.rows) - 1 - self.split) * self.step + 1

    @property
    def scan_y_dist(self) -> float:
        """float: The distance between the rows of the rapid scan pass."""
        return (self.rows[1] - self.rows[0]) / self.step


def plan_rows(
    min_y: float, max_y: float, y_count: int, scan_min_y: float
) -> None | HybridPlan:
    """Plan the Y rows of a hybrid mesh.

    The row distance is floored to the hundredth of a millimeter the same way
    ``bed_mesh`` does, and to a multiple of ``step`` hundredths so the denser
    rapid scan rows land exactly on the mesh rows.

    Args:
        min_y (float): The front edge of the mesh.
        max_y (float): The back edge of the mesh.
        y_count (int): The number of mesh rows.
        scan_min_y (float): The smallest Y the eddy coil can scan.

    Raises:
        BedMeshError: If the eddy coil can not reach at least two mesh rows or
            the rows end up closer than 1 mm.

    Returns:
        None | HybridPlan: None if the eddy coil can reach the whole mesh and
            no tap pass is needed, the plan otherwise.
    """
    if min_y >= scan_min_y:
        return None
    step = 1
    # a larger step makes the row distance a bit smaller, which can push the
    # overlap row one row further back, needing a larger step again, step is
    # at most MIN_SCAN_ROWS - 1 so this settles quickly
    for _ in range(MIN_SCAN_ROWS):
        raw_dist = (max_y - min_y) / (y_count - 1) * 100
        dist = math.floor(round(raw_dist / step, 6)) * step / 100
        if dist / step < 1.0:
            raise BedMeshError("hybrid_bed_mesh: mesh rows are too close together")
        rows = [min_y + i * dist for i in range(y_count)]
        split = next((i for i, y in enumerate(rows) if y >= scan_min_y - 1e-6), None)
        if split is None or split == y_count - 1:
            raise BedMeshError(
                "hybrid_bed_mesh: the eddy probe can reach less than two mesh "
                f"rows (scan_min_y: {scan_min_y:.2f}), use "
                "BED_MESH_CALIBRATE METHOD=tap for this area"
            )
        new_step = max(1, math.ceil((MIN_SCAN_ROWS - 1) / (y_count - 1 - split)))
        if new_step == step:
            return HybridPlan(rows, split, step)
        step = new_step
    raise BedMeshError("hybrid_bed_mesh: could not plan the mesh rows")


def tap_columns(x_count: int, tap_x_count: int) -> list[int]:
    """Pick the mesh columns that are tapped.

    The tapped columns are spread as evenly as possible and always fall on
    mesh columns, so the overlap row can be compared point by point.

    Args:
        x_count (int): The number of mesh columns.
        tap_x_count (int): The number of points to tap per row.

    Returns:
        list[int]: The indices of the tapped columns, ascending.
    """
    tap_x_count = max(2, min(tap_x_count, x_count))
    return sorted(
        {round(j * (x_count - 1) / (tap_x_count - 1)) for j in range(tap_x_count)}
    )


def natural_cubic_spline(
    xs: Sequence[float], ys: Sequence[float], x_new: Sequence[float]
) -> list[float]:
    """Interpolate with a natural cubic spline.

    Args:
        xs (Sequence[float]): The known X values, strictly ascending.
        ys (Sequence[float]): The known Y values.
        x_new (Sequence[float]): The X values to interpolate at, values
            outside of ``xs`` are extrapolated from the end segments.

    Returns:
        list[float]: The interpolated Y values.
    """
    n = len(xs)
    h = [xs[i + 1] - xs[i] for i in range(n - 1)]
    # solve the tridiagonal system for the second derivatives, which are 0 at
    # both ends for a natural spline
    sub = [0.0] * n
    diag = [1.0] * n
    sup = [0.0] * n
    rhs = [0.0] * n
    for i in range(1, n - 1):
        sub[i] = h[i - 1]
        diag[i] = 2 * (h[i - 1] + h[i])
        sup[i] = h[i]
        rhs[i] = 6 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1])
    for i in range(1, n):
        w = sub[i] / diag[i - 1]
        diag[i] -= w * sup[i - 1]
        rhs[i] -= w * rhs[i - 1]
    m = [0.0] * n
    m[-1] = rhs[-1] / diag[-1]
    for i in range(n - 2, -1, -1):
        m[i] = (rhs[i] - sup[i] * m[i + 1]) / diag[i]

    result = []
    for x in x_new:
        i = 0
        while i < n - 2 and x > xs[i + 1]:
            i += 1
        x0, x1 = xs[i], xs[i + 1]
        y0, y1 = ys[i], ys[i + 1]
        hi = h[i]
        result.append(
            m[i] * (x1 - x) ** 3 / (6 * hi)
            + m[i + 1] * (x - x0) ** 3 / (6 * hi)
            + (y0 / hi - m[i] * hi / 6) * (x1 - x)
            + (y1 / hi - m[i + 1] * hi / 6) * (x - x0)
        )
    return result


def merge_rows(
    tap_rows: list[list[float]],
    columns: list[int],
    scan_rows: list[list[float]],
    linear: bool = False,
) -> tuple[list[list[float]], list[float], float]:
    """Merge the tapped front rows with the scanned rows.

    The last tap row and the first scan row are the overlap row. The tap rows
    are corrected by the tap - scan difference on that row, either by its
    average or by a least squares line along X, then interpolated onto every
    mesh column.

    Args:
        tap_rows (list[list[float]]): The tapped Z values, front to back, one
            value per tapped column, the last row is the overlap row.
        columns (list[int]): The mesh column index of each tapped value.
        scan_rows (list[list[float]]): The scanned Z values, front to back,
            one value per mesh column, the first row is the overlap row.
        linear (bool): Fit a line along X to the overlap row differences
            instead of using their average.

    Returns:
        tuple[list[list[float]], list[float], float]: The merged Z matrix
            (front to back, left to right), the correction that was
            subtracted from the tap values at each tapped column and the
            largest remaining tap - scan difference on the overlap row.
    """
    overlap_tap = tap_rows[-1]
    overlap_scan = [scan_rows[0][c] for c in columns]
    diffs = [t - s for t, s in zip(overlap_tap, overlap_scan)]
    offset = sum(diffs) / len(diffs)
    slope = 0.0
    mean_column = sum(columns) / len(columns)
    if linear:
        slope = sum(
            (c - mean_column) * (d - offset) for c, d in zip(columns, diffs)
        ) / sum((c - mean_column) ** 2 for c in columns)
    corrections = [offset + slope * (c - mean_column) for c in columns]
    deviation = max(abs(d - k) for d, k in zip(diffs, corrections))
    x_count = len(scan_rows[0])
    front = [
        natural_cubic_spline(
            columns, [z - k for z, k in zip(row, corrections)], range(x_count)
        )
        for row in tap_rows[:-1]
    ]
    return front + [list(row) for row in scan_rows], corrections, deviation


def check_algorithm(
    algo: str, x_count: int, y_count: int, x_pps: int, y_pps: int
) -> str:
    """Validate the interpolation algorithm the same way ``bed_mesh`` does.

    Args:
        algo (str): The requested algorithm, ``lagrange`` or ``bicubic``.
        x_count (int): The number of mesh columns.
        y_count (int): The number of mesh rows.
        x_pps (int): The interpolated points per segment along X.
        y_pps (int): The interpolated points per segment along Y.

    Raises:
        BedMeshError: If the algorithm can not be used with this mesh.

    Returns:
        str: The algorithm to use.
    """
    if algo not in ("lagrange", "bicubic"):
        raise BedMeshError(f"hybrid_bed_mesh: Unknown algorithm <{algo}>")
    if max(x_pps, y_pps) == 0:
        return "direct"
    max_count = max(x_count, y_count)
    min_count = min(x_count, y_count)
    if algo == "lagrange" and max_count > 6:
        raise BedMeshError(
            "hybrid_bed_mesh: cannot exceed a probe_count of 6 when using "
            f"lagrange interpolation. Probe Count: {x_count}, {y_count}"
        )
    if algo == "bicubic" and min_count < 4:
        if max_count > 6:
            raise BedMeshError(
                "hybrid_bed_mesh: 3 points on one axis with more than 6 on the "
                "other is not permitted with bicubic interpolation. "
                f"Probe Count: {x_count}, {y_count}"
            )
        return "lagrange"
    return algo


def parse_pair(gcmd: GCodeCommand, name: str, default: Sequence, parser: type) -> tuple:
    """Parse a ``X,Y`` pair gcode parameter.

    Args:
        gcmd (GCodeCommand): The gcode command.
        name (str): The parameter name.
        default (Sequence): The value to return if the parameter is missing.
        parser (type): ``float`` or ``int``.

    Raises:
        gcmd.error: If the value is not a pair of numbers.

    Returns:
        tuple: The parsed pair.
    """
    value = gcmd.get(name, None)
    if value is None:
        return tuple(default)
    try:
        pair = tuple(parser(v.strip()) for v in value.split(","))
    except ValueError:
        pair = ()
    if len(pair) != 2:
        raise gcmd.error(f"hybrid_bed_mesh: invalid {name} value '{value}'")
    return pair


class TapPointsHelper(ProbePointsHelper):
    """A ProbePointsHelper that lets the toolhead settle and retries failed taps.

    Klipper starts the probe as soon as the move to the point ends, a tap
    started while the toolhead still rings from that move triggers in mid-air.
    A tap can also fail now and then ("Unable to detect tap: ..."), which would
    otherwise abort the whole mesh, so that point is tapped again.
    """

    settle_time: float = 0.0
    retries: int = 0

    def _move_next(self, probe_num: int) -> None:
        """Move to the next point and wait ``settle_time`` seconds.

        Args:
            probe_num (int): The index of the point to move to.
        """
        super()._move_next(probe_num)
        if self.settle_time:
            self.printer.lookup_object("toolhead").dwell(self.settle_time)

    def start_probe(self, gcmd: GCodeCommand) -> None:
        """Tap every point, retrying a point whose tap couldn't be detected.

        Unlike ``ProbePointsHelper.start_probe`` the results are pulled after
        each point: a failed tap leaves its error in the eddy probe session, so
        a retry needs a new session, and ending a session drops its results.

        Args:
            gcmd (GCodeCommand): The tap command, with ``METHOD=tap``.

        Raises:
            gcmd.error: If a point still fails after ``retries`` retries.
        """
        probe = self.printer.lookup_object("probe")
        self.horizontal_move_z = gcmd.get_float(
            "HORIZONTAL_MOVE_Z", self.default_horizontal_move_z
        )
        self.lift_speed = probe.get_probe_params(gcmd)["lift_speed"]
        self.probe_offsets = probe.get_offsets(gcmd)
        session = probe.start_probe_session(gcmd)
        results = []
        for probe_num, (x, y) in enumerate(self.probe_points):
            for attempt in range(self.retries + 1):
                self._raise_tool(not probe_num and not attempt)
                self._move_next(probe_num)
                try:
                    session.run_probe(gcmd)
                except self.printer.command_error as e:
                    if not str(e).startswith(TAP_ERROR):
                        raise
                    if attempt == self.retries:
                        raise gcmd.error(
                            f"hybrid_bed_mesh: tap at X{x:.2f} Y{y:.2f} failed "
                            f"{attempt + 1} times: {e}"
                        ) from e
                    gcmd.respond_info(
                        f"hybrid_bed_mesh: tap at X{x:.2f} Y{y:.2f} failed: {e}, "
                        f"retrying ({attempt + 1}/{self.retries})"
                    )
                    session.end_probe_session()
                    session = probe.start_probe_session(gcmd)
                else:
                    results.extend(session.pull_probed_results())
                    break
        self._raise_tool()
        session.end_probe_session()
        self.finalize_callback(results)


class HybridBedMesh:
    """The [hybrid_bed_mesh] config section.

    The section name must not start with ``bed_mesh``, bed_mesh loads every
    section with that prefix as a ``[bed_mesh <profile>]`` mesh profile.

    Args:
        config (ConfigWrapper): The config wrapper.
    """

    def __init__(self, config: ConfigWrapper) -> None:
        self.printer: Printer = config.get_printer()
        self.gcode = self.printer.lookup_object("gcode")
        self.bed_mesh: BedMesh = self.printer.load_object(config, "bed_mesh")
        orig_config = self.bed_mesh.bmc.orig_config
        if orig_config["radius"] is not None:
            raise config.error("hybrid_bed_mesh: round beds are not supported")
        self.mesh_min = config.getfloatlist(
            "mesh_min", orig_config["mesh_min"], count=2
        )
        self.mesh_max = config.getfloatlist(
            "mesh_max", orig_config["mesh_max"], count=2
        )
        self.probe_count = config.getintlist(
            "probe_count",
            (orig_config["x_count"], orig_config["y_count"]),
            count=2,
        )
        self.scan_min_y = config.getfloat("scan_min_y", None)
        self.tap_x_count = config.getint("tap_x_count", 7, minval=2)
        self.tap_samples = config.getint("tap_samples", 1, minval=1)
        self.tap_sample_retract_dist = config.getfloat(
            "tap_sample_retract_dist", None, above=0.0
        )
        self.tap_probe_speed = config.getfloat("tap_probe_speed", None, above=0.0)
        self.tap_settle_time = config.getfloat("tap_settle_time", 0.5, minval=0.0)
        self.tap_retries = config.getint("tap_retries", 2, minval=0)
        self.overlap_alignment = config.getchoice(
            "overlap_alignment",
            {name: name for name in OVERLAP_ALIGNMENTS},
            "linear",
        )
        self.max_overlap_offset = config.getfloat("max_overlap_offset", 0.5, above=0.0)
        self.max_overlap_deviation = config.getfloat(
            "max_overlap_deviation", 0.1, above=0.0
        )
        self.position_min_y = config.getsection("stepper_y").getfloat(
            "position_min", 0.0, note_valid=False
        )
        # reads horizontal_move_z and speed from this section, used for the
        # moves between the tapped points
        self.tap_helper = TapPointsHelper(
            config, self._tap_finalize, default_points=[(0.0, 0.0)]
        )
        self.tap_results: list = []
        self.printer.register_event_handler("klippy:connect", self._handle_connect)
        self.gcode.register_command(
            "BED_MESH_HYBRID_CALIBRATE",
            self.cmd_bed_mesh_hybrid_calibrate,
            desc=CALIBRATE_HELP,
        )

    def _handle_connect(self) -> None:
        """Default scan_min_y to the smallest Y the eddy coil can reach."""
        if self.scan_min_y is not None:
            return
        y_offset = self.printer.lookup_object("probe").get_offsets()[1]
        self.scan_min_y = self.position_min_y + max(0.0, y_offset)
        logger.info("hybrid_bed_mesh: scan_min_y defaults to %.3f", self.scan_min_y)

    def _tap_finalize(self, results: list) -> None:
        """Store the tap pass results.

        Args:
            results (list): The ProbeResult of each tapped point.
        """
        self.tap_results = results

    def cmd_bed_mesh_hybrid_calibrate(self, gcmd: GCodeCommand) -> None:
        """Build a bed mesh from an eddy rapid scan and a tapped front strip.

        Args:
            gcmd (GCodeCommand): The gcode command.

        Raises:
            gcmd.error: If any of the passes fail or the overlap row does not
                line up.
        """
        if gcmd.get_int("ADAPTIVE", 0):
            raise gcmd.error(
                "hybrid_bed_mesh: ADAPTIVE=1 is not supported, pass the "
                "adaptive MESH_MIN / MESH_MAX from the slicer instead"
            )
        profile = gcmd.get("PROFILE", "default").strip()
        if not profile:
            raise gcmd.error("Value for parameter 'PROFILE' must be specified")
        min_x, min_y = parse_pair(gcmd, "MESH_MIN", self.mesh_min, float)
        max_x, max_y = parse_pair(gcmd, "MESH_MAX", self.mesh_max, float)
        x_count, y_count = parse_pair(gcmd, "PROBE_COUNT", self.probe_count, int)
        if min(x_count, y_count) < 3:
            raise gcmd.error("hybrid_bed_mesh: PROBE_COUNT must be at least 3")
        if max_x <= min_x or max_y <= min_y:
            raise gcmd.error("hybrid_bed_mesh: invalid min/max points")
        pass_through = {
            name: gcmd.get(name)
            for name in SCAN_PASS_THROUGH
            if gcmd.get(name, None) is not None
        }
        try:
            plan = plan_rows(min_y, max_y, y_count, self.scan_min_y)
        except BedMeshError as e:
            raise gcmd.error(str(e)) from e

        if plan is None:
            gcmd.respond_info(
                "hybrid_bed_mesh: the eddy probe reaches the whole mesh, "
                "skipping the tap pass"
            )
            self._run_scan(
                (min_x, min_y),
                (max_x, max_y),
                (x_count, y_count),
                profile,
                pass_through,
            )
            return

        params, scan_rows = self._scan_pass(
            gcmd, plan, (min_x, max_x), x_count, profile, pass_through
        )
        # the partial mesh is cleared so it does not shift the tap moves
        self.bed_mesh.set_mesh(None)
        columns = tap_columns(x_count, self.tap_x_count)
        tap_rows = self._tap_pass(gcmd, plan, params, columns)
        matrix = self._merge(gcmd, tap_rows, columns, scan_rows)
        algo = gcmd.get("ALGORITHM", self.bed_mesh.bmc.orig_config["algo"])
        try:
            params["algo"] = check_algorithm(
                algo.strip().lower(),
                x_count,
                y_count,
                params["mesh_x_pps"],
                params["mesh_y_pps"],
            )
        except BedMeshError as e:
            raise gcmd.error(str(e)) from e
        params["min_y"] = plan.rows[0]
        params["y_count"] = y_count
        self._set_mesh(gcmd, params, matrix, profile)

    def _scan_pass(
        self,
        gcmd: GCodeCommand,
        plan: HybridPlan,
        x_range: tuple[float, float],
        x_count: int,
        profile: str,
        pass_through: dict[str, str],
    ) -> tuple[dict, list[list[float]]]:
        """Rapid scan the eddy-reachable rows.

        Args:
            gcmd (GCodeCommand): The BED_MESH_HYBRID_CALIBRATE command.
            plan (HybridPlan): The mesh plan.
            x_range (tuple[float, float]): The mesh min and max X.
            x_count (int): The number of mesh columns.
            profile (str): The bed mesh profile name.
            pass_through (dict[str, str]): Extra BED_MESH_CALIBRATE parameters.

        Raises:
            gcmd.error: If the scanned mesh does not match the plan.

        Returns:
            tuple[dict, list[list[float]]]: The scanned mesh parameters and
                the scanned Z values of the mesh rows, front to back.
        """
        # the scan max Y gets a little extra so bed_mesh's flooring of the row
        # distance keeps scan_y_dist
        scan_min_y = plan.rows[plan.split]
        scan_count = plan.scan_row_count
        scan_max_y = scan_min_y + plan.scan_y_dist * (scan_count - 1) + 0.005
        gcmd.respond_info(
            f"hybrid_bed_mesh: 1/2 rapid scanning Y {scan_min_y:.2f} - "
            f"{plan.rows[-1]:.2f} ({scan_count} rows)"
        )
        self._run_scan(
            (x_range[0], scan_min_y),
            (x_range[1], scan_max_y),
            (x_count, scan_count),
            profile,
            pass_through,
        )
        scan_mesh = self.bed_mesh.get_mesh()
        if scan_mesh is None:
            raise gcmd.error("hybrid_bed_mesh: the rapid scan produced no mesh")
        params = dict(scan_mesh.get_mesh_params())
        scan_rows = scan_mesh.get_probed_matrix()[:: plan.step]
        row_count = len(plan.rows) - plan.split
        if (
            len(scan_rows) != row_count
            or len(scan_rows[0]) != x_count
            or abs(params["min_y"] - scan_min_y) > 0.01
            or abs(params["max_y"] - plan.rows[-1]) > 0.01
        ):
            raise gcmd.error(
                "hybrid_bed_mesh: the rapid scan mesh does not match the plan, "
                f"min_y: {params['min_y']:.3f} max_y: {params['max_y']:.3f} "
                f"rows: {len(scan_rows)}, expected min_y: {scan_min_y:.3f} "
                f"max_y: {plan.rows[-1]:.3f} rows: {row_count}"
            )
        return params, scan_rows

    def _tap_pass(
        self,
        gcmd: GCodeCommand,
        plan: HybridPlan,
        params: dict,
        columns: list[int],
    ) -> list[list[float]]:
        """Tap the front rows and the overlap row in a zig-zag.

        Args:
            gcmd (GCodeCommand): The BED_MESH_HYBRID_CALIBRATE command.
            plan (HybridPlan): The mesh plan.
            params (dict): The scanned mesh parameters, for the X positions.
            columns (list[int]): The mesh columns to tap.

        Returns:
            list[list[float]]: The tapped Z values, front to back, left to
                right.
        """
        x_dist = (params["max_x"] - params["min_x"]) / (params["x_count"] - 1)
        xs = [params["min_x"] + c * x_dist for c in columns]
        points = []
        for i in range(plan.split + 1):
            row = [(x, plan.rows[i]) for x in xs]
            points.extend(row if i % 2 == 0 else row[::-1])
        gcmd.respond_info(
            f"hybrid_bed_mesh: 2/2 tapping Y {plan.rows[0]:.2f} - "
            f"{plan.rows[plan.split]:.2f} "
            f"({plan.split + 1} rows x {len(columns)} points)"
        )
        tap_z = self._run_taps(gcmd, points)
        tap_rows = []
        for i in range(plan.split + 1):
            row = tap_z[i * len(columns) : (i + 1) * len(columns)]
            tap_rows.append(row if i % 2 == 0 else row[::-1])
        return tap_rows

    def _merge(
        self,
        gcmd: GCodeCommand,
        tap_rows: list[list[float]],
        columns: list[int],
        scan_rows: list[list[float]],
    ) -> list[list[float]]:
        """Line up the tap rows with the scanned rows and merge them.

        Args:
            gcmd (GCodeCommand): The BED_MESH_HYBRID_CALIBRATE command.
            tap_rows (list[list[float]]): The tapped Z values.
            columns (list[int]): The tapped mesh columns.
            scan_rows (list[list[float]]): The scanned Z values.

        Raises:
            gcmd.error: If OVERLAP_ALIGNMENT is invalid or the tap - scan
                correction exceeds max_overlap_offset.

        Returns:
            list[list[float]]: The merged Z matrix.
        """
        alignment = gcmd.get("OVERLAP_ALIGNMENT", self.overlap_alignment).lower()
        if alignment not in OVERLAP_ALIGNMENTS:
            raise gcmd.error(
                f"hybrid_bed_mesh: invalid OVERLAP_ALIGNMENT '{alignment}', "
                f"use one of {', '.join(OVERLAP_ALIGNMENTS)}"
            )
        matrix, corrections, deviation = merge_rows(
            tap_rows, columns, scan_rows, linear=alignment == "linear"
        )
        offset = sum(corrections) / len(corrections)
        gcmd.respond_info(
            f"hybrid_bed_mesh: overlap row tap - scan offset: {offset:.4f} mm, "
            f"X slope: {corrections[-1] - corrections[0]:.4f} mm across the "
            f"mesh ({alignment}), largest remaining difference: "
            f"{deviation:.4f} mm"
        )
        largest = max(corrections, key=abs)
        if abs(largest) > self.max_overlap_offset:
            raise gcmd.error(
                f"hybrid_bed_mesh: the tap - scan correction of {largest:.4f} mm "
                f"exceeds max_overlap_offset ({self.max_overlap_offset:.4f} "
                "mm), check the tap_threshold, the eddy calibration and the "
                "axis twist compensation"
            )
        if deviation > self.max_overlap_deviation:
            gcmd.respond_info(
                f"hybrid_bed_mesh: WARNING, the tap and scan values on the "
                f"overlap row differ by up to {deviation:.4f} mm after lining "
                f"them up (max_overlap_deviation: "
                f"{self.max_overlap_deviation:.4f} mm)"
            )
        return matrix

    def _set_mesh(
        self,
        gcmd: GCodeCommand,
        params: dict,
        matrix: list[list[float]],
        profile: str,
    ) -> None:
        """Build, activate and save the merged mesh.

        Args:
            gcmd (GCodeCommand): The BED_MESH_HYBRID_CALIBRATE command.
            params (dict): The merged mesh parameters.
            matrix (list[list[float]]): The merged Z matrix.
            profile (str): The bed mesh profile name.

        Raises:
            gcmd.error: If bed_mesh can not build the mesh.
        """
        z_mesh = ZMesh(params, profile)
        try:
            z_mesh.build_mesh(matrix)
        except BedMeshError as e:
            raise gcmd.error(str(e)) from e
        zero_ref = self.bed_mesh.bmc.probe_mgr.get_zero_ref_pos()
        if (
            zero_ref is not None
            and params["min_x"] <= zero_ref[0] <= params["max_x"]
            and params["min_y"] <= zero_ref[1] <= params["max_y"]
        ):
            z_mesh.set_zero_reference(*zero_ref)
        self.bed_mesh.set_mesh(z_mesh)
        gcmd.respond_info("Hybrid Mesh Bed Leveling Complete")
        self.bed_mesh.save_profile(profile)

    def _run_scan(
        self,
        mesh_min: tuple[float, float],
        mesh_max: tuple[float, float],
        probe_count: tuple[int, int],
        profile: str,
        pass_through: dict[str, str],
    ) -> None:
        """Run a BED_MESH_CALIBRATE METHOD=rapid_scan pass.

        Args:
            mesh_min (tuple[float, float]): The mesh min point.
            mesh_max (tuple[float, float]): The mesh max point.
            probe_count (tuple[int, int]): The number of columns and rows.
            profile (str): The bed mesh profile name.
            pass_through (dict[str, str]): Extra BED_MESH_CALIBRATE parameters.
        """
        script = (
            "BED_MESH_CALIBRATE METHOD=rapid_scan "
            f"MESH_MIN={mesh_min[0]:.3f},{mesh_min[1]:.3f} "
            f"MESH_MAX={mesh_max[0]:.3f},{mesh_max[1]:.3f} "
            f"PROBE_COUNT={probe_count[0]},{probe_count[1]} "
            f"PROFILE={profile}"
        )
        for name, value in pass_through.items():
            script += f" {name}={value}"
        self.gcode.run_script_from_command(script)

    def _run_taps(
        self, gcmd: GCodeCommand, points: list[tuple[float, float]]
    ) -> list[float]:
        """Tap the given points with the nozzle.

        Args:
            gcmd (GCodeCommand): The BED_MESH_HYBRID_CALIBRATE command.
            points (list[tuple[float, float]]): The nozzle XY positions.

        Raises:
            gcmd.error: If not every point was tapped.

        Returns:
            list[float]: The bed Z at each point.
        """
        params = {
            "METHOD": "tap",
            "SAMPLES": str(gcmd.get_int("TAP_SAMPLES", self.tap_samples, minval=1)),
        }
        retract_dist = gcmd.get_float(
            "TAP_SAMPLE_RETRACT_DIST", self.tap_sample_retract_dist, above=0.0
        )
        if retract_dist is not None:
            params["SAMPLE_RETRACT_DIST"] = str(retract_dist)
        # METHOD=tap ignores the probe's speed option and defaults to 5mm/s
        probe_speed = gcmd.get_float(
            "TAP_PROBE_SPEED", self.tap_probe_speed, above=0.0
        )
        if probe_speed is not None:
            params["PROBE_SPEED"] = str(probe_speed)
        move_z = gcmd.get_float("TAP_HORIZONTAL_MOVE_Z", None, above=0.0)
        if move_z is not None:
            params["HORIZONTAL_MOVE_Z"] = str(move_z)
        tap_gcmd = self.gcode.create_gcode_command(
            "BED_MESH_HYBRID_CALIBRATE", "BED_MESH_HYBRID_CALIBRATE", params
        )
        self.tap_helper.settle_time = gcmd.get_float(
            "TAP_SETTLE_TIME", self.tap_settle_time, minval=0.0
        )
        self.tap_helper.retries = gcmd.get_int(
            "TAP_RETRIES", self.tap_retries, minval=0
        )
        self.tap_results = []
        self.tap_helper.update_probe_points(points, 1)
        self.tap_helper.start_probe(tap_gcmd)
        if len(self.tap_results) != len(points):
            raise gcmd.error(
                f"hybrid_bed_mesh: tapped {len(self.tap_results)} of "
                f"{len(points)} points"
            )
        return [result.bed_z for result in self.tap_results]


def load_config(config: ConfigWrapper) -> HybridBedMesh:
    """Load the [hybrid_bed_mesh] config section.

    Args:
        config (ConfigWrapper): The config wrapper.

    Returns:
        HybridBedMesh: The hybrid bed mesh instance.
    """
    return HybridBedMesh(config)
