''' Helpers Functions
'''
# Imports
from RootTools.core.Sample import Sample

# Helpers Functions
def singleton(class_):
  instances = {}
  def getinstance(*args, **kwargs):
    if class_ not in instances:
        instances[class_] = class_(*args, **kwargs)
    return instances[class_]
  return getinstance

def merge( pd, totalRunName, listOfRuns, dirs ):
    dirs[ pd + '_' + totalRunName ] = []
    for run in listOfRuns:
        dirs[ pd + '_' + totalRunName ].extend( dirs[ pd + '_' + run ] )

