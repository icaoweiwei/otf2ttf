@echo off
set fontpath=%1

if "%fontpath%"=="" (
  echo Please input dir path!
  pause
  goto end
) 

for %%f in (%fontpath%\*.otf) do (
  echo convert %%f begin
  otfccdump %%f > "%%f.json" && java -jar otf2ttf.jar "%%f.json" | otfccbuild -o "%%f.ttf"
  del /F /Q "%%f.json"
  echo convert %%f end
)

echo Success!

:end