# @Author: Massimo De Mauri <massimo>
# @Date:   2020-12-16T14:06:20+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: load_problem_from_clib.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-17T18:45:41+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function load_problem_from_clib(libPath::String)::Problem

    # load the clib in OpenBB
    id = OpenBB.load_lib(libPath)

    # variables data
    numVars = ccall(dlsym(OpenBB.get_lib_pointer(id),:get_num_vars),Cint,())
    varVals = Vector{Float}(undef,numVars)
    varLoBs = Vector{Float}(undef,numVars)
    varUpBs = Vector{Float}(undef,numVars)
    dscMask = Vector{Cint}(undef,numVars)
    ccall(dlsym(OpenBB.get_lib_pointer(id),:get_vars_data),Cint,(Ptr{Float},Ptr{Float},Ptr{Float},Ptr{Cint}),varVals,varLoBs,varUpBs,dscMask)

    # constraints data
    numCnss = ccall(dlsym(OpenBB.get_lib_pointer(id),:get_num_cnss),Cint,())
    cnsLoBs = Vector{Float}(undef,numCnss)
    cnsUpBs = Vector{Float}(undef,numCnss)
    ccall(dlsym(OpenBB.get_lib_pointer(id),:get_cns_bounds),Cint,(Ptr{Float},Ptr{Float}),cnsLoBs,cnsUpBs)

    # constraints evaluation function
    function evalCnsVal(x::Vector{Float})::Vector{Float}
        out = Vector{Float}(undef,numCnss)
        ccall(OpenBB.dlsym(OpenBB.get_lib_pointer(id),:g_eval),Cint,(Ptr{Float},Ptr{Float}),x,out)
        return out
    end

    # constraints jacobian sparsity
    nnzCJ = ccall(dlsym(OpenBB.get_lib_pointer(id),:get_nnz_Jc),Cint,())

    (rowsCJ,colsCJ) = (Vector{Cint}(undef,nnzCJ),Vector{Cint}(undef,nnzCJ))
    ccall(dlsym(OpenBB.get_lib_pointer(id),:get_sparsity_Jg),Cint,(Ptr{Cint},Ptr{Cint}),rowsCJ,colsCJ)
    (rowsCJ,colsCJ) = (@. Int(rowsCJ+1),@. Int(colsCJ+1))
    cnsJbcSpasity = sparse(rowsCJ,colsCJ,true,numCnss,numVars)

    # evaluate constraints jacobian (in a sparse manner)
    function evalCnsJcb(x::Vector{Float})::SpMatrix{Float}
        values = Vector{Float}(undef,nnzCJ)
        ccall(OpenBB.dlsym(OpenBB.get_lib_pointer(id),:eval_Jg),Cint,(Ptr{Float},Ptr{Float}),x,values)
        return OpenBB.sparse(rowsCJ,colsCJ,values,numCnss,numVars)
    end

    # collect constraints hessians sparsities
    nnzCH = Vector{Cint}(undef,numCnss)
    cnsHesSpasity = Vector{SpMatrix{Bool}}(undef,numCnss)
    names = Vector{String}(undef,numCnss)
    rowsCH = Vector{Union{Vector{Cint},Vector{Int}}}(undef,numCnss)
    colsCH = Vector{Union{Vector{Cint},Vector{Int}}}(undef,numCnss)

    for k in 1:numCnss
        nnzCH[k] = ccall(dlsym(OpenBB.get_lib_pointer(id),Symbol(:get_nnz_Hg,k-1)),Cint,())
        (rowsCH[k],colsCH[k]) = (Vector{Cint}(undef,nnzCH[k]),Vector{Cint}(undef,nnzCH[k]))
        ccall(dlsym(OpenBB.get_lib_pointer(id),Symbol(:get_sparsity_Hg,k-1)),Cint,(Ptr{Cint},Ptr{Cint}),rowsCH[k],colsCH[k])
        if rowsCH[k][1] != -1
            (rowsCH[k],colsCH[k]) = (@. Int(rowsCH[k]+1),@. Int(colsCH[k]+1))
        else
            nnzCH[k]=0
            (rowsCH[k],colsCH[k]) = (Int[],Int[])
        end
        cnsHesSpasity[k] = sparse(rowsCH[k],colsCH[k],true,numVars,numVars)
    end

    # evaluate constraint hessians (sparse)
    function evalCnsHes(x::Vector{Float})::Vector{SpMatrix{Float}}
        out = Vector{SpMatrix{Float}}(undef,numCnss)
        for k in 1:numCnss
            if nnzCH[k] > 0
                values = Vector{Float}(undef,nnzCH[k])
                ccall(OpenBB.dlsym(OpenBB.get_lib_pointer(id),Symbol(:eval_Hg,k-1)),Cint,(Ptr{Float},Ptr{Float}),x,values)
                out[k] = OpenBB.sparse(rowsCH[k],colsCH[k],values,numVars,numVars)
            else
                out[k] = OpenBB.sparse(Int[],Int[],Float[],numVars,numVars)
            end
        end
        return out
    end

    # evaluate objective
    function evalObjVal(x::Vector{Float})::Float
        out = Vector{Float}(undef,1)
        ccall(OpenBB.dlsym(OpenBB.get_lib_pointer(id),:eval_f),Cint,(Ptr{Float},Ptr{Float}),x,out)
        return out[1]
    end

    # objective gradient sparsity
    nnzOG = ccall(dlsym(OpenBB.get_lib_pointer(id),:get_nnz_Gf),Cint,())
    indsOG = Vector{Cint}(undef,nnzOG)
    ~ =  Vector{Cint}(undef,nnzOG)
    ccall(dlsym(OpenBB.get_lib_pointer(id),:get_sparsity_Gf),Cint,(Ptr{Cint},Ptr{Cint}),~,indsOG)
    indsOG = @. Int(indsOG+1)
    objGrdSparsity = sparsevec(indsOG,true,numVars)

    # evaluate objective gradient (sparse)
    function evalObjGrd(x::Vector{Float})::SpVector{Float}
        values = Vector{Float}(undef,nnzOG)
        ccall(OpenBB.dlsym(OpenBB.get_lib_pointer(id),:eval_Gf),Cint,(Ptr{Float},Ptr{Float}),x,values)
        return OpenBB.sparsevec(indsOG,values,numVars)
    end

    # objective hessian sparsity
    nnzOH = ccall(dlsym(OpenBB.get_lib_pointer(id),:get_nnz_Hf),Cint,())
    (rowsOH,colsOH) = (Vector{Cint}(undef,nnzOH),Vector{Cint}(undef,nnzOH))
    ccall(dlsym(OpenBB.get_lib_pointer(id),:get_sparsity_Hf),Cint,(Ptr{Cint},Ptr{Cint}),rowsOH,colsOH)

    if rowsOH[1]!=-1
        (rowsOH,colsOH) = (@. Int(rowsOH+1),@. Int(colsOH+1))
    else
        nnzOH = 0
        (rowsOH,colsOH) = (Int[],Int[])
    end
    objHesSparsity = sparse(rowsOH,colsOH,true,numVars,numVars)

    # evaluate objective hessian (sparse)
    function evalObjHes(x::Vector{Float})::SpMatrix{Float}
        if nnzOH > 0
            values = Vector{Float}(undef,nnzOH)
            ccall(OpenBB.dlsym(OpenBB.get_lib_pointer(id),:eval_Hf),Cint,(Ptr{Float},Ptr{Float}),x,values)
            return OpenBB.sparse(rowsOH,colsOH,values,numVars,numVars)
        else
            return OpenBB.sparse(Int[],Int[],Float[],numVars,numVars)
        end
    end


    # build problem
    if iszero(objHesSparsity) && evalObjVal(zeros(numVars)) == 0.0
        objFun = OpenBB.LinearObjective(L=deepcopy(evalObjGrd(zeros(numVars))))
    else
        objFun = OpenBB.ConvexObjective{SpMatrix{Float},SpVector{Float}}(evalObjVal,evalObjGrd,evalObjHes,objGrdSparsity,objHesSparsity)
    end


    return OpenBB.Problem(varSet=OpenBB.VariableSet(vals=deepcopy(varVals),loBs=deepcopy(varLoBs),upBs=deepcopy(varUpBs),dscIndices=findall(!iszero,dscMask)),
                          cnsSet=OpenBB.ConvexConstraintSet{SpMatrix{Float},SpMatrix{Float}}(evalCnsVal,evalCnsJcb,evalCnsHes,cnsJbcSpasity,cnsHesSpasity,deepcopy(cnsLoBs),deepcopy(cnsUpBs)),
                          objFun=objFun)
end
