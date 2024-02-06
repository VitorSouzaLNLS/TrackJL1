# Accelerator.jl

using Printf
import Base: !=, setproperty!, setfield!, getproperty
using ..Auxiliary
using ..Constants
using ..Elements: Element

export Accelerator, Accelerator!, find_spos, print_accelerator

electron_rest_energy_eV = Constants.electron_rest_energy_eV

mutable struct Accelerator
    energy::Float64
    cavity_state::Auxiliary.BoolState
    radiation_state::Auxiliary.RadiationState
    vchamber_on::Auxiliary.BoolState
    harmonic_number::Int
    lattice::Vector{Element}
    lattice_version::String
    length::Float64
    velocity::Float64
    beta_factor::Float64
    gamma_factor::Float64
end    

function Accelerator!(energy::Real)
    gamma = energy/electron_rest_energy_eV
    beta = sqrt(1 - (1 / gamma^2))
    velocity = beta * light_speed
    accelerator = Accelerator(energy, Auxiliary.off, Auxiliary.radiation_off, Auxiliary.off, 0, Element[], "", 0.0, velocity, beta, gamma)
    return accelerator
end

function isequal(acc1::Accelerator, acc2::Accelerator)
    if (acc1.energy != acc2.energy) return false end
    if (acc1.cavity_state != acc2.cavity_state) return false end
    if (acc1.radiation_state != acc2.radiation_state) return false end
    if (acc1.vchamber_on != acc2.vchamber_on) return false end
    if (acc1.harmonic_number != acc2.harmonic_number) return false end
    if (acc1.lattice_version != acc2.lattice_version) return false end
    #if (acc1.lattice != acc2.lattice) return false end revisar comparador de lattices
    return true
end

function Base.:(!=)(acc1::Accelerator, acc2::Accelerator)
    return !(isequal(acc1, acc2))
end

function update_cavity(accelerator::Accelerator)
    cavity_indices = find_cavity(accelerator)
    for index in cavity_indices
        cav = accelerator.lattice[index]
        if accelerator.cavity_state == on
            cav.properties[:pass_method] = Auxiliary.pm_cavity_pass
        elseif accelerator.cavity_state == off && cav.length == 0.0
            cav.properties[:pass_method] = Auxiliary.pm_identity_pass
        else
            cav.properties[:pass_method] = Auxiliary.pm_drift_pass
        end
    end
end

function setproperty!(accelerator::Accelerator, symbol::Symbol, value)
    if symbol == :energy
        # Custom logic for setting the energy field
        if !(value <= electron_rest_energy_eV)
            #accelerator.energy = value
            setfield!(accelerator, :energy, value)

            gamma = value/electron_rest_energy_eV
            setfield!(accelerator, :gamma_factor, gamma)

            beta = sqrt(1 - (1 / gamma^2))
            setfield!(accelerator, :beta_factor, beta)

            velocity = beta * light_speed
            setfield!(accelerator, :velocity, velocity)
        end
    
    elseif symbol == :cavity_state
        # Custom logic for setting the cavity_state field
        if value == true
            setfield!(accelerator, :cavity_state, on)
        elseif value == false
            setfield!(accelerator, :cavity_state, off)
        elseif value == Auxiliary.on || value == Auxiliary.off
            setfield!(accelerator, :cavity_state, value)
        else
            throw(ArgumentError("Invalid argument for $symbol"))
        update_cavity(accelerator)
        end
    
    elseif symbol == :radiation_state
        # Custom logic for setting the radiation_on field
        setfield!(accelerator, :radiation_state, value)
    
    elseif symbol == :vchamber_on
        # Custom logic for setting the vchamber_on field
        setfield!(accelerator, :vchamber_on, value)
    
    elseif symbol == :harmonic_number
        # Custom logic for setting the harmonic_number field
        setfield!(accelerator, :harmonic_number, value)

    elseif symbol == :lattice
        # Custom logic for setting the lattice field
        setfield!(accelerator, :lattice, value)
        last_index = Int(length(value))+1
        len = find_spos(accelerator, last_index)[1] # closed lattice
        setfield!(accelerator, :length, len)

    elseif symbol == :lattice_version
        # Custom logic for setting the lattice_version field
        setfield!(accelerator, :lattice_version, value)
        @warn("Changing the \"lattice_version\" manually is not recommended..")
    
    elseif symbol == :length # Do nothing -> automatic calculation
        # Custom logic for setting the length field
        # accelerator.length = value
        @warn("Cant manually change the \"length\". Consider changing the accelerator's lattice.")
        # setfield!(accelerator, :length, value)
    
    elseif symbol == :velocity # Do nothing -> automatic calculation
        # Custom logic for setting the velocity field
        #accelerator.velocity = value
        # @warn("Cant manually change the \"velocity\". Consider changing the accelerator's energy.")
        # #update_lorentz_factors(accelerator)

    elseif symbol == :beta_factor # Do nothing -> automatic calculation
        # Custom logic for setting the beta_factor field
        #accelerator.beta_factor = value
        # @warn("Cant manually change the \"beta_factor\". Consider changing the accelerator's energy.")
        # #update_lorentz_factors(accelerator)

    elseif symbol == :gamma_factor # Do nothing -> automatic calculation
        # Custom logic for setting the gamma_factor field
        #accelerator.gamma_factor = value
        # @warn("Cant manually change the \"gamma_factor\". Consider changing the accelerator's energy.")
        #update_lorentz_factors(accelerator)
    else
        throw(ArgumentError("Field $symbol is not a valid field for Accelerator"))
    end
end

function print_accelerator(out::IO, accelerator::Accelerator)
    println(out, "\n--------------- Accelerator ---------------")
    println(out, "energy: $(accelerator.energy) [eV]")
    println(out, "cavity: $(accelerator.cavity_state)")
    println(out, "radiation state: $(accelerator.radiation_state)")
    println(out, "vchamber: $(accelerator.vchamber_on)")
    println(out, "harmonic number: $(accelerator.harmonic_number)")    
    println(out, "length: $(accelerator.length)")
    @printf(out, "velocity: %.8f [m/s]\n", accelerator.velocity)
    @printf(out, "beta factor: %.16f \n", accelerator.beta_factor)
    @printf(out, "gamma factor: %.1f \n", accelerator.gamma_factor)
    if (accelerator.lattice_version != "")
        println(out, "lattice version: $(accelerator.lattice_version)")
    end
    println(out, "-------------------------------------------\n")
end

function find_spos(accelerator::Accelerator, indices::Vector{Int})
    spos = Float64[]
    for index in indices
        if 1 <= index <= length(accelerator.lattice)+1
            spos_elem = sum(element.properties[:length] for element in accelerator.lattice[1:index-1])
            push!(spos, spos_elem)
        else
            error("Invalid index $index for the lattice")
        end
    end
    return spos
end

function find_spos(accelerator::Accelerator, index::Int)
    if 1 <= index <= length(accelerator.lattice)+1
        return sum(element.properties[:length] for element in accelerator.lattice[1:index-1])
    else
        throw(ArgumentError("Invalid index $index for the lattice"))
    end
end

