abstract type AbstractComponent end

abstract type ImplicitComponent <: AbstractComponent end
abstract type ExplicitComponent <: AbstractComponent end

function setName!(comp::AbstractComponent, name)
   comp.name = name
end

function setup!(comp::AbstractComponent)
   return
end

function solveNonlinear!(comp::AbstractComponent, inputs, outputs)
   return
end

function applyNonlinear!(comp::AbstractComponent,
                         inputs,
                         outputs,
                         residuals)
   return
end

"""
   IndepVarComp

Special component that handles the implicit input transformation
"""
mutable struct IndepVarComp{T} <: AbstractComponent
   firstPass::Bool
   targetVals::Dict{String, Base.Vector{T}}
   outputVars::Base.Vector{String}
   outputSizes::Base.Vector{Int}
   name::String

   function IndepVarComp{T}() where {T}
      firstPass = true
      targetVals = Dict{String, Base.Vector{T}}()
      outputVars = Base.Vector{String}()
      outputSizes = Base.Vector{Int}()
      new{T}(firstPass, targetVals, outputVars, outputSizes)
   end
end

function setName!(comp::IndepVarComp, name)
   comp.name = name
end

function copyOutputsToTargets!(comp::IndepVarComp, outputs::Vector)
   for output in comp.outputVars
      comp.targetVals[output] = copy(outputs[output])
   end
   return
end

function solveNonlinear!(comp::IndepVarComp,
                         inputs::Vector,
                         outputs::Vector)
   copyOutputsToTargets!(comp, outputs)
   comp.firstPass = false
   return
end

function applyNonlinear!(comp::IndepVarComp,
                         inputs::Vector,
                         outputs::Vector,
                         residuals::Vector)
   # just do this the first time, in case apply is called before solve
   # otherwise, we assume the target was set by the previous solve
   if comp.firstPass
      copyOutputsToTargets!(comp, outputs)
      comp.firstPass = false
   end
   for output in comp.outputVars
      residuals[output] = outputs[output] - comp.targetVals[output]
   end
   return
end

function solveNonlinear!(comp::ExplicitComponent,
                         inputs::Vector,
                         outputs::Vector)
   compute!(comp, inputs, outputs)
   return
end

function applyNonlinear!(comp::ExplicitComponent,
                         inputs::Vector,
                         outputs::Vector,
                         residuals::Vector)
   # R = compute(inputs, outputs) - outputs
   for output in comp.outputVars
      # subtract output value before compute changes it
      residuals[output] = outputs[output]
   end

   compute!(comp, inputs, outputs)
   for output in comp.outputVars
      # now subtract the newly computed value
      residuals[output] = residuals[output] - outputs[output]
   end
   return
end

function compute!(comp::ExplicitComponent,
                  inputs,
                  outputs)
   return
end
