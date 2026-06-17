# 3J3Q surface ray tracing profile scenario.
# Run Patinae with:
# PATINAE_RT_PROFILE=1 PATINAE_RT_PROFILE_COUNTS=1 PATINAE_RT_PROFILE_STAGES=1 ./target/dev-fast/patinae

load _tests/3J3Q.cif.gz
zoom
as cartoon
set ray_opaque_background, off
set ray_shadow, on
set ray_transparency_shadows, on

# Baseline artifact ray without surface.
ray 1920, 1080

# Surface-heavy path. This should exercise streaming fallback on Metal.
show surface
set surface_transparency, 0.5
ray 1920, 1080

# Repeat at the same camera to expose warm-cache behavior.
ray 1920, 1080
