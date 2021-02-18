# @Author: Massimo De Mauri <massimo>
# @Date:   2021-01-12T16:21:58+01:00
# @Email:  massimo.demauri@protonmail.com
# @Filename: _include.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-12T16:24:39+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# abstract and null types
abstract type SupersolverWorkspace <: AbstractWorkspace end
struct NullSupersolverWorkspace <: SupersolverWorkspace end


# Branch&Bound
include("./branch_and_bound/_include.jl")

# Hybrid-Branch&Bound
include("./hybrid_branch_and_bound/_include.jl")
