using DataStructures

struct Vector{T} <: AbstractVector{T}
   data::Base.Vector{T}
   dataSlices::Dict{String, UnitRange{Int}}
   dataSizes::OrderedDict{String,Int}

   function Vector{T}() where {T}
      data = Base.Vector{T}()
      slices = Dict{String, UnitRange{Int}}()
      sizes = OrderedDict{String, Int}()
      new{T}(data, slices, sizes)
   end
end

Base.IndexStyle(::Type{Vector{T}}) where {T} = IndexLinear()
Base.size(vec::Vector) = size(vec.data)

function Base.setindex!(vec::Vector, value, key::String)
   slice = vec.dataSlices[key]
   vec.data[slice] = value
   return
end

function Base.getindex(vec::Vector, key::String)
   slice = vec.dataSlices[key]
   return @view vec.data[slice]
end

function Base.setindex!(vec::Vector, value, idx::Int)
   vec.data[idx] = value
end

function Base.getindex(vec::Vector, idx::Int)
   return vec.data[idx]
end

function addVar!(vec::Vector, key, size)
   vec.dataSizes[key] = size
   return
end

function allocate!(vec::Vector)
   size = sum(values(vec.dataSizes))
   resize!(vec.data, size)
   # initialize to zero
   vec.data .= zero(eltype(vec.data))

   startIdx = 1
   for (key, size) in vec.dataSizes
      endIdx = startIdx + size-1
      slice = startIdx:endIdx
      vec.dataSlices[key] = slice
      startIdx = endIdx+1
   end

   empty!(vec.dataSizes) # don't need any more
   return
end

function getSlice(vec::Vector, key)
   return vec.dataSlices[key]
end

