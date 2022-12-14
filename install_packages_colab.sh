
#!/bin/bash

mkdir -p /content/drive/MyDrive/installations
cd /content/drive/MyDrive/installations

echo "(1) Install pip packages to /content/drive/MyDrive/installations"
echo "nose mpi4py triangle Pmw pymetis cmocean geopandas fiona pygeos rasterio rasterstats scikit-fmm rtree pyTMD Orinoco"
pip -q install nose mpi4py triangle Pmw pymetis geemap cmocean geopandas fiona pygeos rasterio rasterstats scikit-fmm rtree backports.zoneinfo pyTMD git+https://github.com/simard-landscape-lab/orinoco.git > /dev/null 2>&1

echo "(2) Install gdal"
apt-get -q -y install python-gdal gdal-bin  > /dev/null 2>&1

echo "(3) Install netcdf4"
apt-get -q -y install python-netcdf4  > /dev/null 2>&1

echo "(4) Download anuga_core github repository"
echo "https://github.com/GeoscienceAustralia/anuga_core"
git clone https://github.com/GeoscienceAustralia/anuga_core.git  > /dev/null 2>&1

echo "(5) Install anuga"

cd anuga_core
python setup.py --quiet build  > /dev/null 2>&1
python setup.py --quiet install  > /dev/null 2>&1

cd ../

echo "(7) Completed"
