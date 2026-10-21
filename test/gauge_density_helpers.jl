using Gaugefields
using LinearAlgebra
using Random

const DENSITY_METHOD_KEYS = (
    "plaquette" => "plaquette",
    "clover" => "clover",
    "alexandrou" => "clover improved",
    "bilson_thompson" => "clover improved bilson-thompson",
)

function density_test_field(L; backend=LatticeMatricesBackend(), colors=3, halo=1)
    U = gauge_configuration(L; colors, halo, backend, start=:cold, verbose=0)
    rng = MersenneTwister(210426)
    for site in CartesianIndices(L), μ in 1:4
        A = randn(rng, ComplexF64, colors, colors)
        H = A + A'
        H -= tr(H) / colors * I
        link = exp(0.25im * H)
        for j in 1:colors, i in 1:colors
            U[μ][i, j, Tuple(site)...] = link[i, j]
        end
    end
    Gaugefields.set_wing_U!(U)
    return U
end

function density_test_snapshot(U)
    NC = size(U[1])[1]
    L = size(U[1])[3:end]
    return [U[μ][i, j, Tuple(site)...]
        for μ in 1:4, i in 1:NC, j in 1:NC, site in CartesianIndices(L)]
end

function density_point_reference(U, position, method, m)
    tensor = m.temp_UμνTA
    temps = m._temporary_matrices
    if method == "plaquette"
        return real(QCDMeasurements.calculate_topological_charge_plaq(U, position, tensor, temps))
    end
    q = QCDMeasurements.calculate_topological_charge_clover(U, position, tensor, temps)
    if method == "alexandrou"
        q = QCDMeasurements.calculate_topological_charge_improved(U, position, tensor, q, temps)
    elseif method == "bilson_thompson"
        q = QCDMeasurements.calculate_topological_charge_bilson_thompson(U, position, tensor, q, temps)
    end
    return real(q)
end

function action_point_reference(U, position, method, m)
    if method == "clover"
        nloops = QCDMeasurements.calc_UμνTA!(
            m.temp_UμνTA, position, "clover", U, m._temporary_matrices)
        return -sum(real(tr(m.temp_UμνTA[μ, ν]^2))
                    for μ in 1:4 for ν in (μ + 1):4) / nloops^2
    end
    loops, _ = QCDMeasurements.calc_loopset_μν_name("plaq", 4)
    temps = m._temporary_matrices
    value = 0.0
    for μ in 1:4, ν in (μ + 1):4
        Gaugefields.evaluate_gaugelinks_eachsite!(
            temps[1], loops[μ, ν], U, temps[2:end], position...)
        value += 2 * (size(U[1])[1] - real(tr(temps[1])))
    end
    return value
end
