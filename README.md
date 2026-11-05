# Klipper Hybrid Bed Mesh

A Klipper extra that builds a full bed mesh on printers where the Eddy probe
can't reach the front of the bed.

An Eddy coil mounted behind the nozzle (a positive `y_offset`) runs out of Y
travel before the nozzle does, so a normal `BED_MESH_CALIBRATE` can't measure
the strip at the front of the bed. With a 33.8 mm offset, that's the first
33.8 mm of the bed, often where the first layer starts.

`BED_MESH_HYBRID_CALIBRATE` measures the whole bed in two passes:

1. **Rapid scan:** the Eddy coil rapid scans every row it can reach
   (`BED_MESH_CALIBRATE METHOD=rapid_scan`).
2. **Taps:** the nozzle taps the rows in front of that (`METHOD=tap`). The
   nozzle has no XY offset, so it reaches the front edge. It also taps the
   first scanned row, the **overlap row**.

On the overlap row both passes measure the same points. The tap rows are
shifted by the tap − scan difference there, either a single offset or a
straight line along X, so they line up with the scan. The tap rows are then
interpolated onto every mesh column and merged with the scan into one evenly
spaced mesh. That mesh is saved as a normal `bed_mesh` profile.

```
hybrid_bed_mesh: 1/2 rapid scanning Y 36.41 - 224.87 (19 rows)
hybrid_bed_mesh: 2/2 tapping Y 5.00 - 36.41 (4 rows x 7 points)
hybrid_bed_mesh: overlap row tap - scan offset: 0.0353 mm, X slope: -0.0854 mm across the mesh (linear), largest remaining difference: 0.0319 mm
Hybrid Mesh Bed Leveling Complete
```

## Requirements

- Klipper with Eddy tap support (`PROBE METHOD=tap`, added in January 2026)
  and an Eddy probe with a calibrated `tap_threshold`.
- The Eddy probe behind the nozzle (positive `y_offset`). A probe in front of
  the nozzle or a round bed is not supported.
- A clean nozzle below oozing temperature: filament on the nozzle gives soft
  or early taps.

## Install

```bash
cd ~
git clone https://github.com/eoyilmaz/klipper-hybrid-bed-mesh.git
cd klipper-hybrid-bed-mesh
./install.sh
```

The script links `hybrid_bed_mesh.py` into `~/klipper/klippy/extras/` and
offers to add an `[update_manager]` section to `moonraker.conf`, so Moonraker
can update it. Set `KLIPPER_PATH` if Klipper is somewhere else. It doesn't
change `printer.cfg`: add the `[hybrid_bed_mesh]` section yourself and restart
Klipper.

`./install.sh -u` removes the link and the `[update_manager]` section.

## Configuration

```ini
[hybrid_bed_mesh]
mesh_min: 10, 5
mesh_max: 220, 225
probe_count: 21, 21
scan_min_y: 34
tap_samples: 3
tap_sample_retract_dist: 4
tap_probe_speed: 3
horizontal_move_z: 5
speed: 150
```

[config/example.cfg](config/example.cfg) has the probe, mesh, axis twist and
tap macro sections of a working printer, with comments on why each value was
chosen. [config/start_gcode.gcode](config/start_gcode.gcode) shows the probing
part of a slicer start gcode.

| Option | Default | Description |
|---|---|---|
| `mesh_min` | `[bed_mesh]` `mesh_min` | Front left corner of the whole mesh, in nozzle coordinates. Its Y can be in front of the coil's reach. |
| `mesh_max` | `[bed_mesh]` `mesh_max` | Back right corner of the whole mesh. |
| `probe_count` | `[bed_mesh]` `probe_count` | Columns and rows of the whole mesh, at least 3 each. |
| `scan_min_y` | `stepper_y` `position_min` + probe `y_offset` | The smallest Y the coil scans. Rows in front of it are tapped. |
| `tap_x_count` | 7 | Taps per front row. They're interpolated onto every mesh column with a natural cubic spline. |
| `tap_samples` | 1 | Taps per point, averaged. |
| `tap_sample_retract_dist` | probe default | Lift between the taps of a point. |
| `tap_probe_speed` | 5 mm/s | Tap speed. `METHOD=tap` ignores the probe's `speed` option. |
| `tap_settle_time` | 0.5 | Seconds to wait after each move before tapping. A tap that starts while the toolhead still rings from the move can trigger in mid-air. |
| `tap_retries` | 2 | How many times to tap a point again when its tap can't be detected ("Unable to detect tap: ..."), before the mesh is aborted. |
| `horizontal_move_z` | 5 | Height the taps start from and travel at. |
| `speed` | 50 | Travel speed between the taps, in mm/s. |
| `overlap_alignment` | `linear` | `offset` shifts the tap rows by the average tap − scan difference on the overlap row. `linear` fits a line along X, which also removes an X slope between nozzle and coil. |
| `max_overlap_offset` | 0.5 | Abort if the correction at any tapped column is larger than this, in mm. |
| `max_overlap_deviation` | 0.1 | Warn if tap and scan still differ by more than this after lining them up, in mm. |

The section must not be named `[bed_mesh ...]`: `bed_mesh` treats every such
section as a mesh profile.

## Usage

```
BED_MESH_HYBRID_CALIBRATE [MESH_MIN=x,y] [MESH_MAX=x,y] [PROBE_COUNT=x,y]
    [PROFILE=name] [ALGORITHM=lagrange|bicubic] [MESH_PPS=x,y]
    [HORIZONTAL_MOVE_Z=z] [SCAN_SPEED=v]
    [TAP_SAMPLES=n] [TAP_SAMPLE_RETRACT_DIST=d] [TAP_PROBE_SPEED=v]
    [TAP_SETTLE_TIME=s] [TAP_RETRIES=n] [TAP_HORIZONTAL_MOVE_Z=z]
    [OVERLAP_ALIGNMENT=linear|offset]
```

- `MESH_MIN`, `MESH_MAX`, `PROBE_COUNT` and `PROFILE` work like in
  `BED_MESH_CALIBRATE`, and default to the `[hybrid_bed_mesh]` section.
- `ALGORITHM`, `MESH_PPS`, `HORIZONTAL_MOVE_Z` and `SCAN_SPEED` go to the
  rapid scan.
- The `TAP_*` and `OVERLAP_ALIGNMENT` parameters override the matching
  options for one run.
- `ADAPTIVE=1` isn't supported. Pass the slicer's adaptive mesh area as
  `MESH_MIN` / `MESH_MAX` instead.
- If the whole mesh is within the coil's reach, there's nothing to tap and it
  runs a plain rapid scan.

### Reading the result

The overlap line shows how well the two passes agree:

- **offset:** how much higher the taps read than the scan on the overlap row,
  on average. It's removed from the tap rows.
- **X slope:** how much that difference changes from the left to the right
  end of the mesh. With `linear` it's removed too.
- **largest remaining difference:** what's left after that. A few hundredths
  of a mm is good. Large values mean the taps or the scan are off somewhere;
  check the warning list below.

## Calibration

The merge can only line up the two passes on the overlap row, so both have to
be accurate on their own. In this order:

1. Eddy calibration (`PROBE_EDDY_CURRENT_CALIBRATE`) and the tap threshold
   (`PROBE_EDDY_CURRENT_TAP_CALIBRATE`).
2. Temperature drift calibration (`TEMPERATURE_PROBE_CALIBRATE`). Set
   `resting_z` in `[temperature_probe]` to 4 or more: each sample starts with
   a tap from that height, and the 0.4 default doesn't leave the tap room to
   descend.
3. `Z_TILT_ADJUST` or `QUAD_GANTRY_LEVEL`, if the printer has one.
4. Axis twist compensation along X and Y. Keep `calibrate_end_y` + `y_offset`
   on the bed: the coil has to be over the bed during the tap as well.

## Tap tips

Learned on a BTT Eddy; your printer may differ:

- **Tap at 3 mm/s** (`tap_probe_speed: 3`). At 5 mm/s the results split
  into two separate clusters.
- **Start each tap at least 4 mm above the bed** (`horizontal_move_z: 5`,
  `tap_sample_retract_dist: 4`). Shorter lifts give "no contact found" or
  "invalid depress distance" errors.
- **Clean the nozzle** right before the mesh, and keep it below oozing
  temperature.
- **Heat-soak the bed** before tapping.
- If taps fail in mid-air in one area of the bed only, check the mechanics
  there. On the printer this was written for, front-row taps kept failing
  until a broken Z rod bracket was replaced.

## Development

The planning and merging functions don't need Klipper, and the tests stub the
Klipper modules, so they run without a printer:

```bash
pip install -r requirements-dev.txt
pytest
```

## License

[MIT](LICENSE)
