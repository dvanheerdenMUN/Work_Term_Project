cmake_minimum_required(VERSION 3.5)

project(wt LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#[[
find_package(VTK COMPONENTS
  vtkCommonColor
  vtkCommonCore
  vtkFiltersSources
  vtkFiltersGeometry
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  vtkViewsQt
  vtkIOGeometry
  vtkRenderingAnnotation
  vtkRenderingLabel
  QUIET
)
if (NOT VTK_FOUND)
  message("${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
 ]]#
find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5Charts)

add_executable(wt
  main.cpp
  mainwindow.cpp
  mainwindow.h
  mainwindow.ui
  model.cpp
  model.h
)
target_link_libraries(wt PRIVATE Qt5::Widgets Qt5::Charts Qt5::Core)
#target_link_libraries(wt PRIVATE ${VTK_LIBRARIES})
#target_include_directories(wt PUBLIC /home/s2/Projects-libraries/vtk/Common/Core/)

#target_include_directories(wt PUBLIC "C:/Users/s/Projects/eigen-3.3.7")
target_include_directories(wt PUBLIC "/Users/desiree/Desktop/MITACS/eigen-3.3.7")
target_include_directories(wt PUBLIC "/usr/local/Cellar/boost/1.73.0/include")
