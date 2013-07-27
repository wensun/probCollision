#include <MT/ors.h>
#include <MT/opengl.h>

const char *USAGE=
"\n\
Usage:  ors_edit <ors-filename>\n\
\n\
Iterate between editing the file (with an external editor) and\n\
viewing the model in the OpenGL window (after pressing ENTER).\n\
\n\
Use the number keys 1 2 3 4 5 to toggle display options.\n";

void drawBase(void*){
  glStandardLight();
  glDrawFloor(10,.8,.8,.8);
  glColor(1.,.5,0.);
}

int main(int argn,char **argv){
  cout <<USAGE <<endl;

  const char *file="test.ors";
  if(argn<2){
    cout <<"opening standard file `" <<file <<"'" <<endl;
  }else file=argv[1];

  char *path,*name;
  MT::decomposeFilename(path,name,file);
  chdir(path);
  
  ors::Graph C;
  OpenGL gl;
  gl.add(drawBase,0);
  gl.add(ors::glDrawGraph,&C);
  editConfiguration(name,C,gl);

  return 0;
}
