// 
// Launcher is the main applet class for the rocket flight simulation 
// 
// Launcher uses two other classes:  ControlPanel to handle user input, 
// and Tracker to plot the simulation.  
// 
// Launcher has six methods of its own to handle the central tasks 
// 
//    init()   starts up and initializes the applet, like all applets do 
//    insets() puts a 10-pixel margin on the applet so it looks nice 
// 
//    update() takes input from the control panel and sets simulation variables 
//             this method is activated whenever input is detected from one 
//             of the ControlPanel textFields.  
// 
//    launch() kicks off the simulation, then displays the results through 
//             ControlPanel and Tracker.  
// 
//    trajectory() is the actual routine with the simulation code.  It produces 
//             three arrays with the resulting altitude, velocity and 
//             acceleration for the entire simulation period.  
// 
//    report() takes newly entered time values, updates the displayed  
//             altitude, velocity, and acceleration values, and then replots 
//             the flight curve through Tracker.  
// 


import java.awt.*; 

public class Launcher extends java.applet.Applet 

{  ControlPanel P; 
   Tracker T;  
   double D,Mr,I,Ft,Cd,tc,ChD;   // rocket diameter, rocket mass, impulse, thrust, 
                                 // coef of drag, ejection delay, and chute diameter 
   double[] v,a,y;               // velocity, acceleration, and altitude arrays 

   final double dt = 0.01;       // time increment for simulation 
   final double g  = 9.8;        // acceleration of gravity 
   final double Mc = 0.0107;     // mass of an engine case (Estes standard motors) 
   final int mAS   = 20001;      // size of simulation: simulation time is dt*mAS 
                                 // stands for "max Array Size"  

   public void init() 

   {  int h,w; 
      h = size().height; 
      w = size().width; 

      setLayout(new GridLayout(2,1,10,10)); 

      P = new ControlPanel(this); 
      T = new Tracker(this, h/2-15, w-20); 
          // w-2*inset, h-(2*inset+grid gap)/2

      add(P); 
      add(T); 

      v = new double[mAS]; 
      a = new double[mAS]; 
      y = new double[mAS]; 

   } 

   public Insets insets() 

   {  return new Insets(10,10,10,10); 
   }  

   void update(TextField f) 

   // ControlPanel calls this procedure whenever it detects that an entry 
   // has been made into a textField, or a textField has lost focus 

   {  double entry = 0; 
      boolean goodNumber; 

      goodNumber = f.getText().length() > 0; 
      if (goodNumber && f.isEditable()) 
         try { entry = Double.valueOf(f.getText()).doubleValue(); } 
      catch(NumberFormatException e) 
      {  f.setText(""); 
         System.out.println("Read number error: " + e.getMessage()); 
         goodNumber = false; 
      } 
      if (goodNumber)  
      {  if (f == P.DT) D = entry/1000;       // Get diameter in m 
         else if (f == P.MT) Mr = entry/1000; // Get mass in kg 
         else if (f == P.CdT) Cd = entry;
         else if (f == P.ImT) I  = entry;
         else if (f == P.ThT) Ft = entry;
         else if (f == P.DlT) tc = entry;
         else if (f == P.ChT) ChD = entry/1000; // Get diameter in m 
         if (f.isEditable()) f.setText(String.valueOf(entry)); 
      }
   }  

   void launch() 

   {  int idx, imax; 
      double ymax; 

      trajectory();                  // !!!!! Runs the Simulation !!!!!  

      ymax = 0;                      // Find the max altitude, and 
      imax = 0;                      // the index (= time) where it occurs 
      for (idx=0; idx<mAS; idx++) 
        if (y[idx]>ymax) 
        {  imax = idx; 
           ymax = y[idx];
        }  

      T.update(y,dt,(double)imax/100);        // Plot the altitude trace 

      P.launched = true;                      // Let everyone know that there's 
      P.statusLaunched();                     // valid simulation data available 

      P.TmL.setBackground(Color.black);       // Set up for entry of time values 
      P.TmL.setForeground(Color.yellow); 
      P.TmL.setText("Enter Time"); 
      P.TmT.setEditable(true); 

      P.TmT.setText(String.valueOf((float)Math.round(100*imax*dt)/100) + " s");  
      P.AltT.setText(String.valueOf((float)Math.round(100*y[imax])/100) + " m");  
      P.VelT.setText(String.valueOf((float)Math.round(100*v[imax])/100) + " m/s");  
      P.AccT.setText(String.valueOf((float)Math.round(100*a[imax])/100) + " m/s/s");  
   }  


   void trajectory() 

   // THE HEART OF THE SIMULATION 
   // See web page at http://www.execpc.com/~culp/rockets/rckt_sim.html 
   // for explanation of how this algorithm works 

   {  int idx; 
      double area,A,M,M0,Mp,dM,F,Fd,rho,t,tb,ChA,v0; 

      // First compute constant values 

      Mp = I/800.0;                         // Estimate mass of propellant (black powder) 
      tb = I/Ft;                            // Estimate burn time 
      dM = Mp/tb;                           // Compute mass decrement 
      A  = Math.PI*D*D/4;                   // Compute area of rocket body 
      ChA = Math.PI*ChD*ChD/4;              // Compute area of chute 

      // Set the initial values for loop variables 

      M0 = Mr+Mc+Mp;                                  // initial mass 
      y[0] = v[0] = a[0] = 0;                         // time zero at rest on the ground 

      // Start the simulation loop 

      for(idx=1; idx<mAS; idx++) 
      {  rho = 1.22*Math.pow(0.9,y[idx-1]/1000);  
         v0 = v[idx-1]; 
         t  = idx*dt;  

         M  = t > tb ? Mr + Mc : M0 - dM*t;   
         area = t < tb+tc ? A : Math.max(A,ChA); 
         Fd = v0 == 0 ? 0 : 0.5*rho*Cd*area*v0*v0*v0/Math.abs(v0);  
         F  = t > tb ? -Fd-M*g : Ft-Fd-M*g; 
         a[idx] = F/M; 
         v[idx] = v[idx-1] + a[idx]*dt; 
         y[idx] = y[idx-1] + v[idx]*dt; 
         if (y[idx] < 0) 
         {  y[idx] = 0;  
            v[idx] = 0; 
         }  
      } 
   }

   void report(boolean launched) 

   // ControlPanel calls this procedure whenever it detects that an entry 
   // has been made into the TIME textField 

   {  if (launched) 
      {  double entry = 0;  
         boolean goodNumber; 
         int idx; 

         goodNumber = P.TmT.getText().length() > 0; 
         if (goodNumber) 
            try { entry = Double.valueOf(P.TmT.getText()).doubleValue(); } 
         catch(NumberFormatException e) 
         {  P.TmT.setText(""); 
            System.out.println("Read number error: " + e.getMessage()); 
            goodNumber = false; 
         } 

         if (goodNumber) 
         {  idx = Math.min(Math.round((float)(entry/dt)),mAS-1);  

            P.TmT.setText(String.valueOf((float)Math.round(100*idx*dt)/100) + " s");  
            P.AltT.setText(String.valueOf((float)Math.round(100*y[idx])/100) + " m");  
            P.VelT.setText(String.valueOf((float)Math.round(100*v[idx])/100) + " m/s");  
            P.AccT.setText(String.valueOf((float)Math.round(100*a[idx])/100) + " m/s/s");  
         } 
         T.update(y,dt,Math.min(entry,(double)(mAS-1)/100)); 
       } 
   } 

}  
