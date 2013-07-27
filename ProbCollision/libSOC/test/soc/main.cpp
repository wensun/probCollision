#include <MT/ors.h>
#include <MT/soc.h>
#include <MT/opengl.h>

const char* USAGE="usage: ./x.exe -orsfile test.ors -dynamic 1 -Hcost 1e-3";

int main(int argn,char **argv){
  MT::initCmdLine(argn,argv);
  cout <<USAGE <<endl;

  //-- setup the simulator, swift, opengl, and the SocAbstraction
  ors::Graph ors;
  ors.init(MT::getParameter<MT::String>("orsfile",MT::String("test.ors")));
  //cout <<"read configuration:\n" <<ors <<endl;

  SwiftModule swift;
  swift.init(ors,.5);
  
  OpenGL gl;
  gl.add(glStandardScene);
  gl.add(ors::glDrawGraph,&ors);
  gl.camera.setPosition(5,-10,10);
  gl.camera.focus(0,0,1);
  gl.watch("loaded configuration - press ENTER");

  uint T=200;
  soc::SocSystem_Ors soc;
  if(MT::getParameter<bool>("dynamic",false)){
    soc.initPseudoDynamic(&ors,&swift,&gl,3.,T);
  }else{
    soc.initKinematic(&ors,&swift,&gl,T);
  }
  soc.os=&std::cout;

  //-- setup the control variables (problem definition)
  TaskVariable *pos = new TaskVariable("position",ors, posTVT,"endeff","<t(0 0 .2)>",0,0,ARR());
  pos->setGainsAsAttractor(20,.2);
  pos->y_target = arr(ors.getName("target")->X.p.v,3);
  
  TaskVariable *col = new TaskVariable("collision",ors, collTVT,0,0,0,0,ARR(.05));
  col->setGains(.5,.0);
  col->y_prec=1e-0;
  col->y_target = ARR(0.);

  soc.setTaskVariables(TUPLE(pos,col));
  
  //-- feedback control (kinematic or dynamic) to reach the targets
  arr q,dq,qv;
  soc.getq0(q);
  soc.getqv0(qv);
  for(uint t=0;t<T;t++){
    //soc::bayesianIKControl(soc,dq,0);
    //q += dq;
    if(!soc.dynamic){
      soc::bayesianIKControl2(soc,q,q,0);
      soc.setq(q);
    }else{
      soc::bayesianDynamicControl(soc,qv,qv,0);
      soc.setqv(qv);
    }
    //soc.reportOnState(cout); //->would generate detailed ouput on the state of all variables...
    gl.update(STRING("bayesian Inverse Kinematics: iteration "<<t));
    //gl.watch();
  }
  gl.watch("<press ENTER>");
  
  //-- planning (AICO) to generate an optimal (kinematic) trajectory
  soc.getq0(q);
  soc.setq(q);
  pos->setInterpolatedTargetsEndPrecisions(T,1e-3,1e3,0.,1e3);
  col->setInterpolatedTargetsConstPrecisions(T,1e-2,0.);
  
  q.clear();
  AICO_solver(soc,q,1e-2,.7,.01,0,0);
  ofstream os("z.traj"); q.writeRaw(os); os.close();
  for(;;) soc.displayTrajectory(q,NULL,1,"AICO (planned trajectory)");
  
  return 0;
}
