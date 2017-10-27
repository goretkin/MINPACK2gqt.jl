module FortranHelp
export @aref

@inline function array_pointer_ref(array, indices...)
  pointer(array, sub2ind(array, indices...))
end

macro aref(ex)
  #= tried to make it work without a function invocation
  quote
    arr = $(esc(ex.args[1]))
    # Base.unsafe_convert(Ptr{eltype(arr)}, arr) + sub2ind(arr, idxs...)
    #pointer(arr, sub2ind(arr, esc(ex.args[2:end])...))
    #pointer(arr, $(esc(Expr(:call, sub2ind, quote(arr), ex.args[2:end]...))))
    #pointer(arr, $(Expr(:call, sub2ind, :arr, ex.args[2:end]...)))
  end
  =#
  # copied from Base.@view
  exr = Expr(:call, array_pointer_ref, ex.args...)
  Expr(:&&, true, esc(exr))
end

end
