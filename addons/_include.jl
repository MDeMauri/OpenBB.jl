# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-04T12:29:51+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: addons.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-16T20:52:16+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# check the presence of the MPC addon and load it
if isdir(homeDirectory*"/addons/MPCaddon")
    const WITH_MPC_ADDON = true
    include("./MPCaddon/_include.jl")
else
    const WITH_MPC_ADDON = false
end
