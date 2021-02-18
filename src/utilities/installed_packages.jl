# @Author: Massimo De Mauri <massimo>
# @Date:   2020-12-10T13:30:40+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: installed_packages.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-10T13:30:51+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# checks if a package is installed
function installed()
    deps = dependencies()
    installs = Dict{String, VersionNumber}()
    for (uuid, dep) in deps
        dep.is_direct_dep || continue
        dep.version === nothing && continue
        installs[dep.name] = dep.version
    end
    return installs
end
