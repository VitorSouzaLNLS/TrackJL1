# Track.jl

module Trackjl

export Constants
export Auxiliary
export Elements, Element!, print_element
export Pos
export AcceleratorModule, Accelerator!, find_indices, find_spos, print_accelerator
export Tracking
export read_flatfile!, write_flatfile!
export Sirius

include("src/Pos/posModule.jl")
include("src/Constants/constantsModule.jl")
include("src/Auxiliary/auxiliaryModule.jl")
include("src/Elements/elementsModule.jl")
include("src/Accelerator/acceleratorModule.jl")
include("src/Tracking/trackingModule.jl")
include("src/FlatFile/flatfileModule.jl")
include("src/Models/SI/Sirius.jl")

using .Constants
using .Auxiliary
using .Elements
using .PosModule
using .AcceleratorModule
using .Tracking
using .FlatFile
using .Sirius

end