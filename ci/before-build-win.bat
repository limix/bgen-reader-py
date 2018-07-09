call ci\set-win-path.bat

powershell -Command "(New-Object Net.WebClient).DownloadFile('https://raw.githubusercontent.com/horta/zlib.install/master/install.bat', 'install-zlib.bat')"
powershell -Command "(New-Object Net.WebClient).DownloadFile('https://raw.githubusercontent.com/horta/zstd.install/master/install.bat', 'install-zstd.bat')"
powershell -Command "(New-Object Net.WebClient).DownloadFile('https://raw.githubusercontent.com/horta/almosthere/master/install.bat', 'install-athr.bat')"
powershell -Command "(New-Object Net.WebClient).DownloadFile('https://raw.githubusercontent.com/limix/bgen/master/install.bat', 'install-bgen.bat')"