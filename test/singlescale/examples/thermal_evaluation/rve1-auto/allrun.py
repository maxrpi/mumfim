#!/usr/bin/env python

import os
import sys
import subprocess
import tempfile
import shutil
import re
import numpy as np
from collections import namedtuple

sourcefile_list = [
  "load-env.sh",
  "amsi.yaml",
  "petsc_options",
  "bar.yaml",
  "kraken_auto.xml",
  "prep.sh",
  "run.sh",
]

class gds2command:
  def __init__(self, gdsfile,
              stackfile,
              prefix,
              cellname,
              exe = "/users/bloomm2/gds2model/src/gds2model/bin/gds2model",
              xc=0, yc=0,
              xw=1, yw=1,
              ):
    self.gdsfile = gdsfile
    self.stackfile = stackfile
    self.prefix = prefix
    self.cellname = cellname
    self.xc = xc
    self.yc = yc
    self.xw = xw
    self.yw = yw
    self.exe = exe

  def commandList(self, xc: float, yc: float,
                        xw: float, yw: float,
                        index: int = 0):
    return [
      self.exe,
      "-g", self.gdsfile,
      "-s", self.stackfile,
      "-c", self.cellname,
      "-i",
      "-m", f"{self.prefix}",
      "-l", f"{xc-xw/2}", "-r", f"{xc + xw/2}",
      "-b", f"{yc-yw/2}", "-t", f"{yc + yw/2}",
    ]
    

class meshercommand:
  def __init__(self, 
              prefix,
              exe = "/users/bloomm2/gds2model/src/gds2model/bin/mesher",
              ):
    self.parasolidfile = f"{prefix}.xmt_txt"
    self.simmetrixmodel = f"{prefix}.smd"
    self.simmetrixmesh = f"{prefix}.sms"
    self.mappingfile = f"{prefix}.map"
    self.exe = exe

  def commandList(self):
    return [
      self.exe,
      self.parasolidfile,
      self.simmetrixmodel,
      self.simmetrixmesh,
      self.mappingfile,
    ]
    

class prepcommand:
  def __init__(self, 
              exe = "./prep.sh",
              ):
    self.exe = exe

  def commandList(self):
    return [
      self.exe,
    ]

class mumfimcommand:
  def __init__(self, 
              exe = "./run.sh",
              ):
    self.exe = exe

  def commandList(self):
    return [
      self.exe,
    ]

def populate(dirname):
  for file in sourcefile_list:
    shutil.copy(f"{sourcedirectory}/{file}", f"{dirname}/{file}")
    
def extractKappa(text):
  searchString = "MacroKappa:" + r"\s+([.\-\d]+)"*9
  match = re.search(searchString, text)

  return [match.group(1), match.group(2), match.group(3)], \
        [match.group(4), match.group(5), match.group(6), match.group(7), match.group(8), match.group(9)]

def chunk(n_tasks, n_chunks):
  tasklist = []
  for i in range(n_chunks):
    tasklist.append([j for  j in range(n_tasks * i // n_chunks, n_tasks * (i + 1) // n_chunks)])
  return tasklist

def eprint(msg):
  print(msg, file=sys.stderr)
    
if __name__ == "__main__":

  prefix = "bar"
  cellname = "KRAKEN_CHIP01"
  sourcedirectory = "/lore/bloomm2/mumfim/mumfim/test/singlescale/examples/thermal_evaluation/rve1-auto/"
  rootdir = "/tmp/kappa_working/"

  ## Window to scan over
  xl = 7000.0; xu = 7150.0
  yl = 1850.0; yu = 2000.0

  n_jobs = int(sys.argv[1])
  jobsnum = int(sys.argv[2])
  assert(jobsnum <= n_jobs)

  debugdir = f"{rootdir}/debug"

  os.makedirs(f"{rootdir}", exist_ok=True)
  os.makedirs(f"{debugdir}", exist_ok=True)
  makeModel = gds2command("/users/bloomm2/gds2model/data/kraken.gds",
                          "kraken_auto.xml",
                          f"./{prefix}",
                          cellname=cellname,
                          )  
  makeMesh = meshercommand(prefix=f"./{prefix}")
  prep = prepcommand()
  mumfim = mumfimcommand()

  XW = 2.0; YW = 2.0 # RVE widths
  NX = 17; NY = 16
  dirnames = []
  tasklists = chunk(NX, n_jobs)
  for i in tasklists[jobsnum-1]:
    xc = xl + i * (xu - xl) / (NX - 1)
    for j in range(0, NY):
      yc = yl + j * (yu - yl) / (NY - 1)
      jobindex = i * NX + j
      eprint(f"{i=}, {j=}, {xc=}, {yc=}")

      debug = True
      struct = namedtuple('struct','name')
      tmpdir = tempfile.TemporaryDirectory(prefix=rootdir) if (not debug) else struct(name=f"{debugdir}")
      dirnames.append(tmpdir)
      
      #We want errors to fail, rather than use old files. Only needed if we're not rotating tmpdirs.
      for t_filename in [f"{prefix}.xmt_txt", f"{prefix}.smd", f"{prefix}.smb"]:
        t_filename = f"{tmpdir.name}/{t_filename}"
        if(os.path.exists(t_filename)):
          os.remove(t_filename)
      populate(tmpdir.name)
      
      if debug:
        eprint(f"Running: {' '.join(makeModel.commandList(xc, yc, XW, YW, 0))}")
      output = subprocess.run(makeModel.commandList(xc, yc, XW, YW, 0), cwd=tmpdir.name, stdout=sys.stderr)
      output = subprocess.run(makeMesh.commandList(), cwd=tmpdir.name, stdout=sys.stderr)
      output = subprocess.run(prep.commandList(), cwd=tmpdir.name, stdout=sys.stderr)
      output = subprocess.run(mumfim.commandList(), cwd=tmpdir.name, stdout=subprocess.PIPE)
      

      coord, kappa = extractKappa(output.stdout.decode())

      print(f"{coord=}: {kappa=}")
      
