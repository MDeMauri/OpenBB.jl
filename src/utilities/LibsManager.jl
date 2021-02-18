# @Author: Massimo De Mauri <massimo>
# @Date:   2020-12-02T16:25:54+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: libsManager.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-15T16:12:49+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using Libdl

mutable struct LibsManager
    locations::Vector{String}
    pointers::Vector{Ptr{Nothing}}
    modified::Bool
    function LibsManager()::LibsManager
        return new(String[],Ptr{Nothing}[],false)
    end
end

function list_libs()::Vector{String}
    return libsManager.locations
end

function get_lib_pointer(id::Int)::Ptr{Nothing}
    return libsManager.pointers[id]
end


function load_lib(location::String)::Int
    push!(libsManager.locations,location)
    if location != ""
        push!(libsManager.pointers,dlopen(location))
        try ccall(dlsym(libsManager.pointers[end],:init),Cint,()) catch end
    else
        push!(libsManager.pointers,C_NULL)
    end
    libsManager.modified=true
    return length(libsManager.locations)
end

function load_libs(locations::Vector{String})::Vector{Int}
    len_ = length(libsManager.locations)
    @. load_lib(locations)
    return collect(len_+1:len_+length(locations))
end


function unload_lib(id::Int)::Nothing
    dlclose(libsManager.pointers[id])
    libsManager.pointers[id]=C_NULL
    libsManager.locations[id]=""
    libsManager.modified=true
    return
end

function unload_libs(ids::AbstractArray)::Nothing
    @. unload_lib(ids)
    return
end

function unload_all_libs()::Nothing
    for k in 1:length(libsManager.locations)
        if isempty(libsManager.locations[k])
            dlclose(libsManager.pointers[k])
            libsManager.pointers[k] = C_NULL
        end
        libsManager.modified=true
    end
    return
end

function remove_all_libs()::Nothing
    unload_all_libs()
    empty!(libsManager.locations)
    empty!(libsManager.pointers)
    return
end


# define an object used to bring clibs to remote processes
libsManager = LibsManager()
