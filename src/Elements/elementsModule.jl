module Elements

include("Elements.jl")

export marker, corrector, hcorrector, vcorrector, drift, matrix, rbend, dipole, quadrupole, sextupole, rfcavity, kickmap

# Define functions for Element
function marker(fam_name::String)
    element = Element!(fam_name)
    element.properties[:pass_method] = pm_identity_pass
    return element
end

function corrector(fam_name::String, hkick::Float64, vkick::Float64)
    element = Element!(fam_name)
    element.properties[:pass_method] = pm_corrector_pass
    element.properties[:hkick] = hkick
    element.properties[:vkick] = vkick
    return element
end

function hcorrector(fam_name::String, hkick::Float64)
    element = Element!(fam_name)
    element.properties[:pass_method] = pm_corrector_pass
    element.properties[:hkick] = hkick
    element.properties[:vkick] = 0.0
    return element
end

function vcorrector(fam_name::String, vkick::Float64)
    element = Element!(fam_name)
    element.properties[:pass_method] = pm_corrector_pass
    element.properties[:hkick] = 0.0
    element.properties[:vkick] = vkick
    return element
end

function drift(fam_name::String, length::Float64)
    element = Element!(fam_name)
    element.properties[:pass_method] = pm_drift_pass
    element.properties[:length] = length
    return element
end

function matrix(fam_name::String, length::Float64)
    element = Element!(fam_name)
    element.properties[:pass_method] = pm_matrix_pass
    element.properties[:length] = length
    return element
end

function rbend(fam_name::String, angle::Float64, angle_in::Float64, angle_out::Float64,
                gap::Float64, fint_in::Float64, fint_out::Float64, polynom_a::Vector{Float64},
                polynom_b::Vector{Float64}, K::Float64, S::Float64, nr_steps::Int)
    element = Element!(fam_name)
    element.properties[:pass_method] = pm_bnd_mpole_symplectic4_pass
    element.properties[:angle] = angle
    element.properties[:angle_in] = angle_in
    element.properties[:angle_out] = angle_out
    element.properties[:gap] = gap
    element.properties[:fint_in] = fint_in
    element.properties[:fint_out] = fint_out
    element.properties[:polynom_a] = polynom_a
    element.properties[:polynom_b] = polynom_b
    element.properties[:polynom_b][2] = K
    element.properties[:polynom_b][3] = S
    element.properties[:nr_steps] = nr_steps
    return element
end

function dipole(fam_name::String, angle::Float64, angle_in::Float64, angle_out::Float64,
                           gap::Float64, fint_in::Float64, fint_out::Float64, polynom_a::Vector{Float64},
                           polynom_b::Vector{Float64}, K::Float64, S::Float64, nr_steps::Int)
    return rbend(fam_name, angle, angle_in, angle_out, gap, fint_in, fint_out, polynom_a, polynom_b, K, S, nr_steps)
end

function quadrupole(fam_name::String, K::Float64, length::Float64, nr_steps::Int)
    element = Element!(fam_name)
    element.properties[:pass_method] = pm_str_mpole_symplectic4_pass
    element.properties[:polynom_a] = [0.0, 0.0, 0.0] # copy(default_polynom)
    element.properties[:polynom_b] = [0.0, 0.0, 0.0] # copy(default_polynom) 
    element.properties[:polynom_b][2] = K
    element.properties[:nr_steps] = nr_steps
    element.properties[:length] = length
    return element
end

function sextupole(fam_name::String, S::Float64, nr_steps::Int)
    element = Element!(fam_name)
    element.properties[:pass_method] = pm_str_mpole_symplectic4_pass
    element.properties[:polynom_a] = [0.0, 0.0, 0.0] # copy(default_polynom)
    element.properties[:polynom_b] = [0.0, 0.0, 0.0] # copy(default_polynom) 
    element.properties[:polynom_b][3] = S
    element.properties[:nr_steps] = nr_steps
    return element
end

function rfcavity(fam_name::String, frequency::Float64, voltage::Float64, phase_lag::Float64)  
    element = Element!(fam_name)
    element.properties[:pass_method] = pm_cavity_pass
    element.properties[:frequency] = frequency
    element.properties[:voltage] = voltage
    element.properties[:phase_lag] = phase_lag
    return element
end

# Revisar implementacao
function kickmap(fam_name::String, kicktable_idx::Int, nr_steps::Int, rescale_kicks::Float64)
    element = Element!(fam_name)
    element.properties[:pass_method] = pm_kickmap_pass
    element.properties[:nr_steps] = nr_steps
    element.properties[:kicktable_idx] = kicktable_idx
    element.properties[:rescale_kicks] = rescale_kicks
    return element
end


end