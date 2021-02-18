# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-12T15:39:53+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: heuristics.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-17T14:44:55+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

####################### Functions For No Heuristics Case #################################
# builds a null heuristic that cannot be ever called
function setup(problem::Problem,settings::NullBBheuristicsSettings,subsolverSettings::AbstractSettings)::NullBBheuristicsWorkspace
    return NullBBheuristicsWorkspace()
end

# fake outdated marking operation
function make_outdated!(heuristicsWS::NullBBheuristicsWorkspace)::Nothing
    return
end



include("./triggering_conditions.jl")
include("./simple_rounding_heuristics.jl")
include("./feasibility_pump.jl")
