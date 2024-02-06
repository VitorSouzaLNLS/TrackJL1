# Elements.jl

using ..Auxiliary

export Element, Element!, print_element

# Define Element type
mutable struct Element
    fam_name::String
    properties::Dict{Symbol, Any}
end

# Define default constructor for Element
"""aaa"""
function Element!(fam_name::String)
    return Element(fam_name, Dict(:pass_method => pm_identity_pass, :length => 0.0, :vchamber => vchamber_rectangle))
end

# Other functions
function print_element(out::IO, element::Element)
    custom_order = [:pass_method, :nr_steps, :angle, :angle_in, :angle_out, :gap, :fint_in, :fint_out, :polynom_a, :polynom_b, :frequency, :voltage, :phase_lag, :kicktable_idx, :rescale_kicks, :hkick, :vkick]
    
    keys_sorted = [key in keys(element.properties) ? key : nothing for key in custom_order]
    keys_sorted = filter(x -> x !== nothing, keys_sorted)

    name = element.fam_name
    println(out, "$name")
    for key in keys_sorted
        value = element.properties[key]
        println(out, "$key : $value")
    end
end
