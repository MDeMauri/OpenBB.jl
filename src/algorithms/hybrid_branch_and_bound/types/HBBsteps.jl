# @Author: Massimo De Mauri <massimo>
# @Date:   2020-12-10T12:19:38+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: HBBsteps.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-10T12:20:45+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# abstract types
abstract type HBBnlpStepWorkspace <: AbstractWorkspace end
abstract type HBBnlpStepSettings <: AbstractSettings end
abstract type HBBmipStepWorkspace <: AbstractWorkspace end
abstract type HBBmipStepSettings <: AbstractSettings end
