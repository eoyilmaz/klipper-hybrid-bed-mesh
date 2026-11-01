; The probing part of a slicer start gcode (OrcaSlicer placeholders), heat the
; bed before this and keep the nozzle below oozing temperature until the mesh
; is done.

G28 ; home all axis

Z_TILT_ADJUST
G28 Z

WIPE_NOZZLE                          ; your nozzle cleaning macro, a dirty nozzle gives soft or early taps
G0 Z10 F1200                         ; lift clear of the bed before travelling
G0 X117.5 Y117.5 F{travel_speed*60}  ; the tap position, same as [bed_mesh] zero_reference_position
TAP_Z_CALIBRATE                      ; corrects Klipper's Z with the average of 3 taps

; rapid scans with the eddy and taps the rows in front of its reach
BED_MESH_HYBRID_CALIBRATE MESH_MIN={adaptive_bed_mesh_min[0]},{adaptive_bed_mesh_min[1]} MESH_MAX={adaptive_bed_mesh_max[0]},{adaptive_bed_mesh_max[1]} ALGORITHM=[bed_mesh_algo] PROBE_COUNT={bed_mesh_probe_count[0]},{bed_mesh_probe_count[1]}
