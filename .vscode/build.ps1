param($file, $fileDirname, $fileBasenameNoExtension)

$p_inc = python -c "import sysconfig; print(sysconfig.get_path('include'), end='')"
$n_inc = python -c "import numpy; print(numpy.get_include(), end='')"
$p_lib = python -c "import sys, os; print(os.path.join(sys.prefix, 'libs'), end='')"
$p_ver = python -c "import sys; print('python' + str(sys.version_info.major) + str(sys.version_info.minor), end='')"

& "C:\msys64\mingw64\bin\g++.exe" -fdiagnostics-color=always -g -std=c++20 "$file" -I"$p_inc" -I"$n_inc" -L"$p_lib" -l"$p_ver" -o "$fileDirname\$fileBasenameNoExtension.exe"