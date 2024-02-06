# Track.jl

module Trackjl

export Constants
export Auxiliary
export Element, Element!, print_element
export Elements
export Pos
export AcceleratorModule, Accelerator, Accelerator!
export Tracking
export read_flatfile!, write_flatfile!

include("src/Pos/posModule.jl")
include("src/Constants/constantsModule.jl")
include("src/Auxiliary/auxiliaryModule.jl")
include("src/Elements/elementsModule.jl")
include("src/Accelerator/acceleratorModule.jl")
include("src/Tracking/trackingModule.jl")
include("src/FlatFile/flatfileModule.jl")

using .Constants
using .Auxiliary
using .Elements
using .AcceleratorModule
using .PosModule:Pos
using .AcceleratorModule: Accelerator, Accelerator!
using .Tracking
using .FlatFile: read_flatfile!, write_flatfile!

end