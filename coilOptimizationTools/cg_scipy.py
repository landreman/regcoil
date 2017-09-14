from evaluateRegcoil import coilFourier, evaluateFunctionRegcoil, evaluateGradientRegcoil
import numpy as np
import sys
from regcoilScan import readVariable, namelistLineContains
import os
from scipy import optimize
import scipy

# nmax and mmax are maximum Fourier modes to read from nescin file
def cg_scipy(regcoil_input_file):

  # Create coilFourier object
  nmax_sensitivity = readVariable("nmax_sensitivity","int",regcoil_input_file,required=True)
  mmax_sensitivity = readVariable("mmax_sensitivity","int",regcoil_input_file,required=True)
  nescin_filename = readVariable("nescin_filename","string",regcoil_input_file,required=True)
  # Order to use for norm of the gradient
  norm = readVariable("norm","int",regcoil_input_file,required=True)
  # Tolerance for norm of the gradient
  gtol = readVariable("gtol","float",regcoil_input_file,required=True)
  # Number of iterations
  maxiter = readVariable("maxiter","int",regcoil_input_file,required=True)
  # Max n in nescin file
  nmax = readVariable("nmax","int",regcoil_input_file,required=True)
  # max m in nescin file
  mmax = readVariable("mmax","int",regcoil_input_file,required=True)
  
  f = open(nescin_filename,'r')
  nescinObject = coilFourier(nmax,mmax,regcoil_input_file)
  nescinObject.set_Fourier_from_nescin(nescin_filename)

  x0 = nescinObject.omegas_sensitivity

  f = evaluateFunctionRegcoil
  fprime = evaluateGradientRegcoil
  args = (nescinObject,)
  full_output = True
  disp = True
  retall = True
  [xopt,fopt,func_calls,grad_calls,warnflag,allvecs] = scipy.optimize.fmin_cg(f=f,x0=x0,fprime=fprime,args=args,gtol=gtol,norm=norm,maxiter=maxiter,full_output=full_output,disp=disp,retall=retall)

if __name__ == "__main__":
  cg_scipy(sys.argv[1])

