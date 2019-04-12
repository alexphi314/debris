%%  Alex Philpott
%% Connect to STK and return a new TLE
function [line1,line2] = tle_from_stk(r,v,num)
%TLE_FROM_STK Given state vector in ECI, return TLE
% [line1,line2] = tle_from_stk(r,v,num)

%% Connect to STK
uiApplication = actxGetRunningServer('STK11.application');
root = uiApplication.Personality2;
root.NewScenario('Debris_TLE_Scenario');

%% Initialize Satellite
satellite = root.CurrentScenario.Children.New('eSatellite','debris');
%satellite.SetPropagatorType('ePropagatorSGP4');
propagator = satellite.Propagator;
%propagator.PropagationFrame = 'eCoordinateSystemJ2000';
propagator.InitialState.Representation.AssignCartesian('eCoordinateSystemJ2000',r(1),r(2),r(3),v(1),v(2),v(3));
propagator.Propagate;

%% Generate TLE
today = datevec(datetime('now'));
tom = datevec(datetime('tomorrow'));

today = [today(1:3),17,00,00];
tom = [tom(1:3),15,00,00];
tod_str = datestr(today,'dd mmm yyyy HH:MM:SS.FFF');
tom_str = datestr(tom,'dd mmm yyyy HH:MM:SS.FFF');
cmd = sprintf('GenerateTLE */Satellite/debris Sampling "%s" "%s" 300.0 "%s" %i 21 0.0003 SGP4',...
              tod_str,tom_str,tod_str,num);
          
root.ExecuteCommand([cmd]);
DP = satellite.DataProviders.Item('Element Set').Exec;
cell = DP.DataSets.GetDataSetByName('Two Line').GetValues;

line1 = string(cell(1));
line2 = string(cell(2));
root.CloseScenario;
end

