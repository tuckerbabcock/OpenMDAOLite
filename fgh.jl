include("OpenMDAOLite.jl")
import .OpenMDAOLite 
om = OpenMDAOLite

mutable struct F <: om.ExplicitComponent
   inputVars::Base.Vector{String}
   inputSizes::Base.Vector{Int}
   outputVars::Base.Vector{String}
   outputSizes::Base.Vector{Int}
   name::String

   function F()
      inputVars = Base.Vector{String}()
      inputSizes = Base.Vector{Int}()
      outputVars = Base.Vector{String}()
      outputSizes = Base.Vector{Int}()
      new(inputVars, inputSizes, outputVars, outputSizes)
   end
end

function om.setup!(comp::F)
   comp.inputVars = ["x", "h"]
   comp.inputSizes = [2, 1]
   comp.outputVars = ["f"]
   comp.outputSizes = [1]
   return
end

function om.compute!(comp::F, inputs, outputs)
   outputs["f"] = sum(inputs["x"].^2) .+ exp.(inputs["h"])
   return
end

mutable struct R_g <: om.ImplicitComponent
   inputVars::Base.Vector{String}
   inputSizes::Base.Vector{Int}
   outputVars::Base.Vector{String}
   outputSizes::Base.Vector{Int}
   name::String

   function R_g()
      inputVars = Base.Vector{String}()
      inputSizes = Base.Vector{Int}()
      outputVars = Base.Vector{String}()
      outputSizes = Base.Vector{Int}()
      new(inputVars, inputSizes, outputVars, outputSizes)
   end
end

function om.setup!(comp::R_g)
   comp.inputVars = ["x", "f"]
   comp.inputSizes = [2, 1]
   comp.outputVars = ["y_g"]
   comp.outputSizes = [1]
   return
end

function om.solveNonlinear!(comp::R_g, inputs, outputs)
   outputs["y_g"] = sqrt.((prod(inputs["x"]) .+ inputs["f"]).^2)
   # println("solve out ", out)
   # println("out y_g: ", outputs["y_g"], " in x:", inputs["x"], " in f: ", inputs["f"])
   # println()

   return
end

function om.applyNonlinear!(comp::R_g,
                            inputs,
                            outputs,
                            residuals)
   residuals["y_g"] = outputs["y_g"].^2 .- (prod(inputs["x"]) .+ inputs["f"]).^2
   # println("apply res ", res)
   # println("out y_g: ", outputs["y_g"], " in x:", inputs["x"], " in f: ", inputs["f"])
   # println()

   return
end

mutable struct H <: om.ExplicitComponent
   inputVars::Base.Vector{String}
   inputSizes::Base.Vector{Int}
   outputVars::Base.Vector{String}
   outputSizes::Base.Vector{Int}
   name::String

   function H()
      inputVars = Base.Vector{String}()
      inputSizes = Base.Vector{Int}()
      outputVars = Base.Vector{String}()
      outputSizes = Base.Vector{Int}()
      new(inputVars, inputSizes, outputVars, outputSizes)
   end
end

function om.setup!(comp::H)
   comp.inputVars = ["x", "y_g"]
   comp.inputSizes = [2, 1]
   comp.outputVars = ["h"]
   comp.outputSizes = [1]
   return
end

function om.compute!(comp::H, inputs, outputs)
   outputs["h"] = 1.0 ./ (prod(inputs["x"]) .* inputs["y_g"])
   return
end

type = Complex{Float64}
# type = Float64
prob = om.Problem{type}()

indeps = om.addSubsystem!(prob, "indeps", om.IndepVarComp{type}())
indeps.outputVars = ["x"]
indeps.outputSizes = [2.0]

om.addSubsystem!(prob, "f_comp", F())
om.addSubsystem!(prob, "R_g_comp", R_g())
om.addSubsystem!(prob, "h_comp", H())

# Note: all variables are promoted to the top level group and connections are made by variable name
om.setup!(prob)

# print(prob.outputVec)

# initial conditions
prob["x"] = [1.0, 2.0]
prob["h"] = [3.0]

# # simple nonlinear block Gauss-Seidel solver (fixed point iteration)
for i in 0:4
   om.runSolveNonlinear!(prob)
   om.runApplyNonlinear!(prob)

   println("iteration $i residuals")
   for key in keys(prob.resVec.dataSlices)
      println(key, " ", prob.resVec[key])
   end
   println("###########################")
   println()
end

println("Final solution")
for key in keys(prob.outputVec.dataSlices)
    println(key, " ", prob.outputVec[key])
end

println()
println()
##################
# UDE Derivatives
##################

# allocate memory for ∂R/∂y
vecLen = length(prob.resVec)
J_partial = zeros(vecLen, vecLen)

# use complex-step to compute columns of J
###########################################
oVecData = prob.outputVec.data
rVecData = prob.resVec.data
for i in 1:vecLen 
   @. oVecData[:] = complex(real(oVecData), 0.0) # this allocates memory :(
   oVecData[i] += complex(0.0, 1e-50)

   om.runApplyNonlinear!(prob)
   J_partial[:,i] = imag(rVecData) ./ 1e-50
end

# use the vector slices to figure out where to put the seed
###########################################################
xSlice = om.getSlice(prob.outputVec, "x")
seed = zeros(vecLen)
seed[xSlice[1]] = 1

# forward mode UDE
##################
dall_dx0 = J_partial \ seed

# in reverse mode, you put the seed into the output you are interested in 
#########################################################################
fSlice = om.getSlice(prob.outputVec, "f")
seed = zeros(vecLen)
seed[fSlice[1]] = 1
# reverse mode UDE
df_dall = J_partial' \ seed # one row of J_total

# CS total check:
###########################################################
# reconverge the system
J_total_CS = zeros(vecLen, vecLen)
for i in 1:vecLen
   @. oVecData[:] = complex(real(oVecData), 0.0) # this allocates memory :(
   oVecData[i] += complex(0.0, 1e-50)
    
   # reconverge the model for every complex-step
   for j in 1:5
      om.runSolveNonlinear!(prob)
   end

   J_total_CS[:,i] = imag(oVecData) / 1e-50
end

println("dall_dx0", dall_dx0) # one column of J_total
println("dall_dx0 cs check", J_total_CS[:,xSlice[1]])
println()
println("df_dall", df_dall) # one row of J_total
println("df_dall cs check", J_total_CS[fSlice[1],:])
println()