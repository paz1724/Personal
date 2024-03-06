function tap_table = Get_tap_table(channel)

% returns a channel's tap table in the format:
% tap table: [ tap_time_in_ns   tap_power_in_dB ]
%  TS 36.101  tables B.2.1

if strcmpi(channel,'EVEH_A') || strcmpi(channel,'EVA')
    tap_table = [0 30 150 310 370 710 1090 1730 2510
        0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9]';
end
if strcmpi(channel,'EPED_A') || strcmpi(channel,'EPA')
    tap_table = [0 30 70 90 110 190 410
        0 -1 -2 -3 -8 -17.2 -20.8]';
end
if strcmpi(channel,'NEW')
    tap_table = [...
        0,  7000,   14000,  2000;...
        0, 7,     -3,      -2]';
end
if strcmpi(channel,'ETU')
    tap_table = [...
        0,  50, 120,    200,    230,    500,    1600,   2300,   5000,   50000;...
        -1, -1, -1,     0,      0,      0,     -3,     -5,      -7,     7]';
end
if strcmpi(channel,'ETU_like')
    tap_table = [0 50 120 200 230 500 1600 2300 5000
        -1 -1 -1   0      0      0     -3     -5      -17 ]';
end
if strcmpi(channel,'ETU_like1')
    tap_table = [0 50 120 200 230 500 1600 5000
        -1 -1 -1   0      0      0     -3     -17 ]';
end
if strcmpi(channel,'EDS')
    tap_table = [  0          30         150         310         370        1090       12490       12520       12640       12800       12860       13580       27490       27520  27640       27800       27860       28580 ;
        0   -1.5000   -1.4000   -3.6000   -0.6000   -7.0000  -10.0000  -11.5000  -11.4000  -13.6000  -10.6000  -17.0000  -20.0000  -21.5000  -21.4000  -23.6000  -20.6000 -27 ]';
end
if strcmpi(channel,'EDS_2_clusters')
    tap_table = [  0          30         150         310         370        1090       12490       12520       12640       12800       12860       13580       ;
        0   -1.5000   -1.4000   -3.6000   -0.6000   -7.0000  -10.0000  -11.5000  -11.4000  -13.6000  -10.6000  -17.0000  ]';
end
if strcmpi(channel,'EDS_Full')
    slope = -10/13000;
    tap_table(:,1) = [0:500:27000]';
    tap_table(:,2) = slope*tap_table(:,1);
    %         tap_pwr = slope*taps
    %         tap_table = [  0          30         150         310         370        1090       12490       12520       12640       12800       12860       13580       27490       27520  27640       27800       27860       28580 ;
    %                                     0   -1.5000   -1.4000   -3.6000   -0.6000   -7.0000  -10.0000  -11.5000  -11.4000  -13.6000  -10.6000  -17.0000  -20.0000  -21.5000  -21.4000  -23.6000  -20.6000 -27 ]';
end
if strcmpi(channel,'FLAT')|| strcmpi(channel,'AWGN')||strcmpi(channel,'ONE_MINUS_ONE')||strcmpi(channel,'EYE')||strcmpi(channel,'ONE_MINUS_J')
    tap_table = [0  0];
    tap_table(:,1) = tap_table(:,1);
end
if strcmpi(channel,'test')
    tap_table(:,1) = [0:9]*100;
    tap_table(:,2) = 0;
end
%Channel for CQI tests (36.101@B.2.4)
if strcmpi(channel,'TEST_CQI')
    tap_table = [0 450 ; 0 0 ]';
end

if strcmpi(channel,'Repeter_channel_13012020')
    load('H:\yossip\SVN\Lab\tools\CE\Repeter_channel_13012020')
    tap_table=[t;taps]';
end

if strcmpi(channel,'two_taps')
tap_table = [0 17e3 ; 0 0 ]';
end
        