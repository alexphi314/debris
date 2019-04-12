%%  Alex Philpott
%% Connect to STK and return a new TLE
function [line1,line2] = tle_from_stk(r,v,num,date)
%TLE_FROM_STK Given state vector in ECI, return TLE
% [line1,line2] = tle_from_stk(r,v,num,date)

%% Connect to STK
uiApplication = actxGetRunningServer('STK11.application');
root = uiApplication.Personality2;
scenario = root.Children.New('eScenario','Debris_TLE_Scenario');
date_str = datestr(date,'dd mmm yyyy HH:MM:SS.FFF');
scenario.SetTimePeriod(date_str,'+24hr');

%% Initialize Satellite
satellite = root.CurrentScenario.Children.New('eSatellite','debris');
%satellite.SetPropagatorType('ePropagatorSGP4');
propagator = satellite.Propagator;
%propagator.PropagationFrame = 'eCoordinateSystemJ2000';
propagator.InitialState.Representation.AssignCartesian('eCoordinateSystemJ2000',r(1),r(2),r(3),v(1),v(2),v(3));
propagator.Propagate;

%% Generate TLE
startt = date;
startt(4) = startt(4)+1;
endt = [startt(1:2),startt(3)+1,date(4)-1,00,00];
tod_str = datestr(startt,'dd mmm yyyy HH:MM:SS.FFF');
tom_str = datestr(endt,'dd mmm yyyy HH:MM:SS.FFF');
cmd = sprintf('GenerateTLE */Satellite/debris Sampling "%s" "%s" 300.0 "%s" %i 21 0.0003 SGP4',...
              tod_str,tom_str,date_str,num);
          
root.ExecuteCommand([cmd]);
DP = satellite.DataProviders.Item('Element Set').Exec;
cell = DP.DataSets.GetDataSetByName('Two Line').GetValues;

line1 = string(cell(1));
line2 = string(cell(2));
root.CloseScenario;
end

