# Workaround: ThermoFunConfig.cmake includes its targets before calling
# find_package(autodiff/ChemicalFun), so the imported targets are missing.
# Loading them first (via -DCMAKE_PROJECT_INCLUDE=...) fixes find_package(ThermoFun).
find_package(ChemicalFun CONFIG QUIET)
find_package(autodiff CONFIG QUIET)
find_package(spdlog CONFIG QUIET)
