classdef LoudSpeaker
    
    properties
        numOfMics       % Number of microphones in the array
        centerPos       % Center of the array [x y]
        micPos          % Position of the microphones [x1 y1; ... ; [xM yM]]
        micSignals      % Microphone Data [S1; ... ; SM]
        radius          % Radius of the array
        dim             % Dimension of the array
        micRIRs         % Struct for the microphone RIRs
        c               % Speed of sound
        fs              % Sampling Frequency
        dist2Src        % Range towards the source
    end
    
    methods
        function obj = LoudSpeaker(centerPos, numOfMics, radius, dim)
            % CLASS CONSTRUCTOR: Init the LoudSpeaker object containing a UCA.
            % function obj = LoudSpeaker(centerPos, numOfMics, radius, dim)
            % --IN--
            % centerPos      : Position of the center of the array
            % numOfMics      : Number of microphones in the array
            % radius         : Radius of the circular array
            % dim (optional) : Dimension of the simulation (2D or 3D [Default])
            % --OUT--
            % obj : Loudspeaker object
            if nargin < 4
                obj.dim = 3;
            else
                if 3 < dim || 2 > dim 
                    error('Dimension should be greater than 2 and smaller than 3')
                else
                    obj.dim = dim;
                end
            end
            obj.radius = radius;
            obj.centerPos = centerPos;
            obj.numOfMics = numOfMics;
            obj.micPos = obj.genMicPos();
            obj.c = 344.8;                 % Hardcoded speed of sound [m/s]
        end
        
        function micPos = genMicPos(obj)
            % METHOD: Generates the position of the microphones in the UCA
            %  micPos = genMicPos(obj)
            % --OUT--
            % micPos : Position of the microphones [x1 y1;....;xM yM]
            phi = 2*pi/obj.numOfMics * (0:obj.numOfMics - 1)';
            pos = obj.radius * [cos(phi) sin(phi)...
                zeros(obj.numOfMics,1)];
            micPos = repmat(obj.centerPos,obj.numOfMics,1) + ...
                pos;
        end
        
        function plotArray(obj)
            % METHOD: Plots the Array
            % plotArray(obj)
            plot(obj.micPos(:,1), obj.micPos(:,2), 'or')
            hold on, 
            plot(obj.centerPos(1), obj.centerPos(2), 'x')
        end
        
        function obj = genRIRs(obj,srcPos,roomDim,T60,fs)
            % METHOD : Generates the RIRs for each of the microphones in
            % the array
            % obj = genRIRs(obj,srcPos,roomDim,T60,fs)
            % --IN--
            % srcPos    :  Position of the source
            % roomDim   : Dimensions of the room [x y z]
            % T60       : Reberveration time in sec [s] for T60
            % fs        : Sampling frequency
            % --OUT--
            % obj : Updated object
            obj.fs = fs;
            nsample = fs*T60;            % Number of samples
            mtype = 'omnidirectional';   % Type of microphone
            order = -1;                  % ?1 equals maximum reflection order!
            orientation = [0 0];         % Microphone orientation [azimuth elevation] in radians
            hp_filter = 1;               % Enable high-pass filter

            for k = 1:obj.numOfMics
                obj.micRIRs(k,:) =...
                    rir_generator(obj.c, obj.fs, obj.micPos(k,:), srcPos, roomDim, T60,...
                    nsample, mtype, order, obj.dim, orientation, hp_filter);
            end
        end
        
        function obj = genMicSignals(obj, srcSignal)
            % METHOD : Generates the signals received at the microphone in
            % the array
            % obj = genMicSignals(obj, srcSignal)
            % --IN--
            % srcSignal : Temporal signal transmitted by the source
            % --OUT--
            % obj : Updated object
            for k = 1:obj.numOfMics
                obj.micSignals(k,:) = filter(obj.micRIRs(k,:),1,srcSignal);
            end
        end
        
        function obj = genNewSim(obj, srcSignal, srcPos, roomDim, T60, fs)
            % METHOD : Generates a new simulation creating the RIRs for
            % each microphone in the UCA and the appropiate received
            % signal.
            % obj = genNewSim(obj, srcSignal)
            % --IN--
            % srcSignal : Temporal signal transmitted by the source
            % srcPos    : Position of the source
            % roomDim   : Dimensions of the room [x y z]
            % T60       : Reberveration time in sec [s] for T60
            % fs        : Sampling frequency
            % --OUT--
            % obj : Updated object
            obj = obj.genRIRs(srcPos, roomDim, T60, fs);
            obj = obj.genMicSignals(srcSignal);
            obj.dist2Src = norms(...
                obj.micPos - repmat(srcPos,obj.numOfMics,1),2,2);
        end
    end
    
end