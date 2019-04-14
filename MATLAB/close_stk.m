%%  Alex Philpott
function [] = close_stk()
uiApplication = actxGetRunningServer('STK11.application');
root = uiApplication.Personality2;
root.CloseScenario;
end
