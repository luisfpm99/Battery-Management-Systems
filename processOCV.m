function model = processOCV(data,cellID)
% processOCV.m
% Reconstructed study version based on Gregory L. Plett Lesson 2.2.5.
% The lecture slides show code highlights, not the complete original file.
% Therefore, some variables and processing steps used by the original
% implementation are not defined in the slides.
%
% Purpose:
% Estimate the static OCV-SOC-temperature relationship of a battery cell
% from laboratory test data.

filetemps = [data.temp];                     % Extract all experiment temperatures.
filetemps = filetemps(:);                    % Convert temperatures into a column vector.
numtemps = length(filetemps);                % Count the number of temperature datasets.

ind25 = find(filetemps == 25);               % Find the dataset measured at 25 degC.

if isempty(ind25)                            % Check whether 25 degC data exists.
    error('Must have a test at 25 degC');    % Stop if no 25 degC test is available.
end                                          % End of if block.

not25 = find(filetemps ~= 25);               % Find datasets measured at temperatures other than 25 degC.
data25 = data(ind25);                        % Store the 25 degC dataset separately.


% -------------------------------------------------------------------------
% COULOMBIC EFFICIENCY AND CAPACITY
% -------------------------------------------------------------------------

totDisAh = data25.script1.disAh(end) + ...   % Final discharged Ah from script 1.
           data25.script2.disAh(end) + ...   % Add discharged Ah from script 2.
           data25.script3.disAh(end) + ...   % Add discharged Ah from script 3.
           data25.script4.disAh(end);        % Add discharged Ah from script 4.

totChgAh = data25.script1.chgAh(end) + ...   % Final charged Ah from script 1.
           data25.script2.chgAh(end) + ...   % Add charged Ah from script 2.
           data25.script3.chgAh(end) + ...   % Add charged Ah from script 3.
           data25.script4.chgAh(end);        % Add charged Ah from script 4.

eta25 = totDisAh/totChgAh;                  % Compute coulombic efficiency at 25 degC.

data25.script1.chgAh = ...                   % Correct charged Ah in script 1.
    data25.script1.chgAh*eta25;

data25.script2.chgAh = ...                   % Correct charged Ah in script 2.
    data25.script2.chgAh*eta25;

data25.script3.chgAh = ...                   % Correct charged Ah in script 3.
    data25.script3.chgAh*eta25;

data25.script4.chgAh = ...                   % Correct charged Ah in script 4.
    data25.script4.chgAh*eta25;

Q25 = data25.script1.disAh(end) + ...        % Begin capacity estimate at 25 degC.
      data25.script2.disAh(end) - ...        % Add discharged Ah from script 2.
      data25.script1.chgAh(end) - ...        % Subtract corrected charged Ah from script 1.
      data25.script2.chgAh(end);             % Subtract corrected charged Ah from script 2.


% -------------------------------------------------------------------------
% ESTIMATE VOLTAGE CHANGES ASSOCIATED WITH R0
% -------------------------------------------------------------------------

indD = find(data25.script1.step == 2);       % Find samples belonging to the slow-discharge step.

IR1Da = data25.script1.voltage(indD(1)-1) - ... % Voltage immediately before discharge starts.
        data25.script1.voltage(indD(1));         % Voltage immediately after discharge starts.

IR2Da = data25.script1.voltage(indD(end)+1) - ... % Voltage immediately after discharge ends.
        data25.script1.voltage(indD(end));         % Voltage immediately before discharge ends.

indC = find(data25.script3.step == 2);       % Find samples belonging to the slow-charge step.

IR1Ca = data25.script3.voltage(indC(1)) - ... % Voltage immediately after charge starts.
        data25.script3.voltage(indC(1)-1);     % Voltage immediately before charge starts.

IR2Ca = data25.script3.voltage(indC(end)) - ... % Voltage immediately before charge ends.
        data25.script3.voltage(indC(end)+1);    % Voltage immediately after charge ends.

IR1D = min(IR1Da,2*IR2Ca);                  % Limit discharge correction at the beginning.
IR2D = min(IR2Da,2*IR1Ca);                  % Limit discharge correction at the end.

IR1C = min(IR1Ca,2*IR2Da);                  % Limit charge correction at the beginning.
IR2C = min(IR2Ca,2*IR1Da);                  % Limit charge correction at the end.


% -------------------------------------------------------------------------
% ADJUST DISCHARGE VOLTAGE CURVE
% -------------------------------------------------------------------------

blend = (0:length(indD)-1)/(length(indD)-1); % Create values that progress from 0 to 1.

IRblend = IR1D + (IR2D-IR1D)*blend(:);      % Linearly vary the IR correction across discharge.

disV = data(k).script1.voltage(indD) + ...  % Measured discharge voltage.
       IRblend;                             % Add IR correction to approximate OCV.

disZ = 1 - data25.script1.disAh(indD)/Q25;  % Estimate SOC during discharge.

disZ = disZ + (1-disZ(1));                  % Force the discharge curve to begin at 100% SOC.


% -------------------------------------------------------------------------
% ADJUST CHARGE VOLTAGE CURVE
% -------------------------------------------------------------------------

blend = (0:length(indC)-1)/(length(indC)-1); % Create values that progress from 0 to 1.

IRblend = IR1C + (IR2C-IR1C)*blend(:);      % Linearly vary the IR correction across charge.

chgV = data25.script3.voltage(indC) - ...   % Measured charge voltage.
       IRblend;                             % Subtract IR correction to approximate OCV.

chgZ = data25.script3.chgAh(indC)/Q25;      % Estimate SOC during charge.

chgZ = chgZ - chgZ(1);                     % Force the charge curve to begin at 0% SOC.


% -------------------------------------------------------------------------
% COMPENSATE FOR STEADY-STATE RESISTANCE
% -------------------------------------------------------------------------

deltaV50 = interp1(chgZ,chgV,0.5) - ...     % Interpolated charge voltage at 50% SOC.
           interp1(disZ,disV,0.5);          % Minus interpolated discharge voltage at 50% SOC.

ind = find(chgZ < 0.5);                    % Select the lower half of the charge SOC range.

vChg = chgV(ind) - ...                     % Charge voltage values.
       chgZ(ind)*deltaV50;                  % Correct them using the 50% SOC voltage difference.

zChg = chgZ(ind);                          % Store the corresponding charge SOC values.

ind = find(disZ > 0.5);                    % Select the upper half of the discharge SOC range.

vDis = flipud(disV(ind) + ...              % Correct discharge voltage.
       (1-disZ(ind))*deltaV50);             % Apply the 50% SOC voltage correction.

zDis = flipud(disZ(ind));                  % Reverse discharge SOC order.

rawocv = interp1([zChg;zDis], ...           % Combine charge and discharge SOC values.
                 [vChg;vDis], ...           % Combine corrected charge/discharge voltages.
                 SOC, ...                   % Evaluate OCV at the desired SOC points.
                 'linear', ...              % Use linear interpolation.
                 'extrap');                 % Allow extrapolation outside measured points.

filedata(ind25).rawocv = rawocv;           % Store the approximate OCV curve at 25 degC.
filedata(ind25).temp = data25.temp;         % Store the corresponding temperature.


% -------------------------------------------------------------------------
% COMBINE DATA FROM ALL TEMPERATURES
% -------------------------------------------------------------------------

Vraw = [];                                 % Empty array for OCV curves.
temps = [];                                % Empty array for temperatures.

for k = 1:numtemps                         % Loop through all temperature datasets.

    if filedata(k).temp > 0                % Use datasets satisfying this condition.

        Vraw = [Vraw; ...                  % Append one OCV curve.
                filedata(k).rawocv];

        temps = [temps; ...                % Append its corresponding temperature.
                 filedata(k).temp];

    end                                    % End of if block.

end                                        % End of for loop.


% -------------------------------------------------------------------------
% FINAL LEAST-SQUARES FIT
% -------------------------------------------------------------------------

X = [ones(size(temps)),temps] \ Vraw;      % Least-squares fit of OCV versus temperature.

model.OCV0 = X(1,:);                       % Store the base OCV relationship.
model.OCVrel = X(2,:);                     % Store the temperature-dependent OCV term.
model.SOC = SOC;                           % Store the SOC vector used by the model.

end                                        % End processOCV function.
