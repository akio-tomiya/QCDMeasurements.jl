# Standalone CPU/GPU density regression, runnable with one or more MPI ranks.
# Select the JACC backend through the active environment's preferences.
using JACC
JACC.@init_backend
using MPI, QCDMeasurements, Test
MPI.Initialized() || MPI.Init()
include("gauge_density_backends.jl")
