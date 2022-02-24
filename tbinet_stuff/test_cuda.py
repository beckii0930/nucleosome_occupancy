import numpy as np
import scipy.io
from sklearn import metrics
import pandas as pd
import os
os.environ['THEANO_FLAGS'] = "device=cuda0,force_device=True,floatX=float32"
from aesara_theano_fallback import aesara as theano
#import theano

print(theano.config.device)
