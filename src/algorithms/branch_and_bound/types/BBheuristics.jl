# @Author: Massimo De Mauri <massimo>
# @Date:   2020-04-27T18:56:37+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBheuristics.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-04T14:25:14+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# Abstract types
abstract type BBheuristicsWorkspace <: AbstractWorkspace end
abstract type BBheuristicsSettings <: AbstractSettings end

# null types
struct NullBBheuristicsWorkspace <: BBheuristicsWorkspace end
struct NullBBheuristicsSettings <: BBheuristicsSettings end
