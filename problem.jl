struct Problem{T}
   components::Base.Vector{AbstractComponent}
   outputVec::Vector{T}
   scratchVec::Vector{T}
   resVec::Vector{T}

   function Problem{T}() where {T}
      components = Base.Vector{AbstractComponent}()
      outputVec = Vector{T}()
      scratchVec = Vector{T}()
      resVec = Vector{T}()
      new{T}(components, outputVec, scratchVec, resVec)
   end
end

function Base.setindex!(prob::Problem, value, key::String)
   prob.outputVec[key] = value
   return
end

function Base.getindex(prob::Problem, key::String)
   return prob.outputVec[key]
end

function setup!(prob::Problem{T}) where {T}
   for comp in prob.components
      # println(comp.name)
      setup!(comp)
      # println(comp.outputVars)

      for (var, size) in zip(comp.outputVars, comp.outputSizes)
         addVar!(prob.outputVec, var, size)
         addVar!(prob.scratchVec, var, size)
         addVar!(prob.resVec, var, size)
      end
   end
   allocate!(prob.outputVec)
   allocate!(prob.scratchVec)
   allocate!(prob.resVec)
   return
end

function addSubsystem!(prob::Problem, name, comp)
   setName!(comp, name)
   push!(prob.components, comp)
   return comp
end
      
function runSolveNonlinear!(prob::Problem)
   for comp in prob.components
      solveNonlinear!(comp, prob.outputVec, prob.outputVec)
   end
   return
end

function runApplyNonlinear!(prob::Problem)
   prob.scratchVec.data[:] = prob.outputVec.data
   for comp in prob.components
      applyNonlinear!(comp, prob.scratchVec, prob.outputVec, prob.resVec)
   end
   return
end