# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-23T15:59:54+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: HBBheuristics.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-21T10:58:59+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# Abstract types
abstract type HBBheuristicsWorkspace <: AbstractWorkspace end
abstract type HBBheuristicsSettings <: AbstractSettings end

# null types
struct NullHBBheuristicsWorkspace <: HBBheuristicsWorkspace outdated::Bool end
struct NullHBBheuristicsSettings <: HBBheuristicsSettings end
