# *vesselpy* - Hydrostatic solver

[1] Clone the repository:
```
cd $HOME && mkdir VESSELPY && cd VESSELPY
git clone https://github.com/francescosalvadore/vesselpy.git
```

[2] Install requirements (virtual env is recommended):
```
python3 -m venv pyenv
. pyenv/bin/activate
cd vesselpy
pip3 install -r requirements.txt
```

[3] Add vesselpy bin to PATH environment variable (to `.bashrc`):
```
export PATH=$HOME/VESSELPY/vesselpy/bin:$PATH
```

[4] Download and install ParaView:
```
apt install paraview
```

[5] Set ParaView folder in `load_env.sh`, e.g, e.g..:
```
export PARAVIEW_FOLDER="/opt/paraviewopenfoam510/lib/python3.10/site-packages"
```

[6] Run the example:
```
cd .. && mkdir run && cd run
cp -r ../vesselpy/examples/w3/ .
cd w3
vesselpy
```

[7] Check the results:
```
cat output/1/results.dat
cat output/1/stab.dat
```

