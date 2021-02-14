clc
clear all

addpath mfcc

Tw = 25;           % analysis frame duration (ms)
Ts = 10;           % analysis frame shift (ms)
alpha = 0.97;      % preemphasis coefficient
R = [ 300 3700 ];  % frequency range to consider
M = 20;            % number of filterbank channels
C = 12;            % number of cepstral coefficients
L = 22;            % cepstral sine lifter parameter


% N = length(x);
% t = N/fs;

% hamming window (see Eq. (5.2) on p.73 of [1])
hamming = @(N)(0.54-0.46*cos(2*pi*[0:N-1].'/(N-1)));

s = 10^20;


[a,fs1]= audioread('Test_word.m4a');
x = resample(a,16000,fs1);


AudioInputFolder='C:\Users\lvars\OneDrive\Documents\ASU\fall 19\speech analysis\Assignment 3\training';
WORDS = dir(fullfile(AudioInputFolder,'*.m4a'));

counter = 0;
length(WORDS)


for j=1:length(WORDS)
    AudioFile=fullfile(AudioInputFolder,WORDS(j).name);
    [b,fs2] = audioread(AudioFile);                      % read audio and return the sampled data and its sampling rate
    y = resample(b,16000,fs2,j);

    counter = counter +1

    if length(x) > length(y)
        p = length(x) - length(y);
        B = padarray(y,p,'post');
        A = x;
    elseif length(y) > length(x)
        q = length(y) - length(x);
        A = padarray(x,q,'post');
        B = y;
    else
        A = x;
        B = y;
    end

    size(A)
    size(B)

    % Feature extraction (feature vectors as columns)
[ MFCCs1, FBEs1, frames1,eframes1 ] = ...
              mfcc( A, fs1, Tw, Ts, alpha, hamming, R, M, C, L );

% % Plot cepstrum over time
% figure('Position', [30 100 800 200], 'PaperPositionMode', 'auto', ...
%      'color', 'w', 'PaperOrientation', 'landscape', 'Visible', 'on' );
%
% imagesc( [1:size(MFCCs1,2)], [0:C-1], MFCCs1 );
% axis( 'xy' );
% xlabel( 'Frame index' );
% ylabel( 'Cepstrum index' );
% title( 'Mel frequency cepstrum' );

for i = 1:size(MFCCs1,2)
    if and(i > 1,i<size(MFCCs1,2))
        MFCCd1(:,i) = (MFCCs1(:,i+1) - MFCCs1(:,i-1))/2;
        eframesd1(i) = (eframes1(i+1) - eframes1(i-1))/2;
    else
        MFCCd1(:,i) = MFCCs1(:,i);
        eframesd1(i) = eframes1(i);
    end



end

for k = 1:size(MFCCs1,2)
    if and(k > 1,k<size(MFCCs1,2))
         MFCCdd1(:,k) = (MFCCd1(:,k+1) - MFCCd1(:,k-1))/2;
         eframesdd1(k) = (eframesd1(k+1) - eframesd1(k-1))/2;
    else
        MFCCdd1(:,k) = MFCCd1(:,k);
        eframesdd1(k) = eframesd1(k);
    end
end

 if length(MFCCs1) > length(MFCCd1)
        g = length(MFCCs1) - length(MFCCd1);
        MFCCd = padarray(MFCCd1,[0,g],'post');
        MFCCs = MFCCs1;
    elseif length(MFCCd1) > length(MFCCs1)
        e = length(MFCCd1) - length(MFCCs1);
        MFCCs = padarray(MFCCs1,[0,e],'post');
        MFCCd = MFCCd1;
    else
        MFCCs = MFCCs1;
        MFCCd = MFCCd1;
 end

 if length(eframes1) > length(eframesd1)
        r = length(eframes1) - length(eframesd1);
        eframesd = padarray(eframesd1,[0,r],'post');
        eframes = eframes1;
    elseif length(eframesd1) > length(eframes1)
        t = length(eframesd1) - length(eframes1);
        eframes = padarray(eframes1,[0,t],'post');
        eframesd = eframesd1;
    else
        eframes = eframes1;
        eframesd = eframesd1;
 end

MFCCfeats1 = [MFCCs;eframes;MFCCd;MFCCdd1;eframesd;eframesdd1];
%Converting NaN to zeros
MFCCfeats1(isnan(MFCCfeats1))=0;
%Converting the matrix into a vector for distance calculation
MFCCfeats1_vector = MFCCfeats1(:);

% Feature extraction (feature vectors as columns)
[ MFCCs2, FBEs2, frames2,eframes2 ] = ...
              mfcc( y, fs2, Tw, Ts, alpha, hamming, R, M, C, L );

% % Plot cepstrum over time
% figure('Position', [30 100 800 200], 'PaperPositionMode', 'auto', ...
%      'color', 'w', 'PaperOrientation', 'landscape', 'Visible', 'on' );
%
% imagesc( [1:size(MFCCs2,2)], [0:C-1], MFCCs2 );
% axis( 'xy' );
% xlabel( 'Frame index' );
% ylabel( 'Cepstrum index' );
% title( 'Mel frequency cepstrum' );

for i = 1:size(MFCCs2,2)
    if and(i > 1,i<size(MFCCs2,2))
        MFCCd2(:,i) = (MFCCs2(:,i+1) - MFCCs2(:,i-1))/2;
        eframesd2(i) = (eframes2(i+1) - eframes2(i-1))/2;
    else
        MFCCd2(:,i) = MFCCs2(:,i);
        eframesd2(i) = eframes2(i);
    end
end

for i = 1:size(MFCCs2,2)
    if and(i > 1,i<size(MFCCs2,2))
         MFCCdd2(:,i) = (MFCCd2(:,i+1) - MFCCd2(:,i-1))/2;
         eframesdd2(i) = (eframesd2(i+1) - eframesd2(i-1))/2;
    else
        MFCCdd2(:,i) = MFCCd2(:,i);
        eframesdd2(i) = eframesd2(i);
    end
end

 if length(MFCCs2) > length(MFCCd2)
        h = length(MFCCs2) - length(MFCCd2);
        MFCCdn = padarray(MFCCd2,[0,h],'post');
        MFCCsn = MFCCs2;
    elseif length(MFCCd2) > length(MFCCs2)
        l = length(MFCCd2) - length(MFCCs2);
        MFCCsn = padarray(MFCCs2,[0,l],'post');
        MFCCdn = MFCCd2;
    else
        MFCCsn = MFCCs2;
        MFCCdn = MFCCd2;
 end

 if length(eframes2) > length(eframesd2)
        v = length(eframes2) - length(eframesd2);
        eframesdn = padarray(eframesd2,[0,v],'post');
        eframesn = eframes2;
    elseif length(eframesd2) > length(eframes2)
        w = length(eframesd2) - length(eframes2);
        eframesn = padarray(eframes2,[0,w],'post');
        eframesdn = eframesd2;
    else
        eframesn = eframes2;
        eframesdn = eframesd2;
 end

% Combining all the features in a matrix
MFCCfeats2 = [MFCCsn;eframesn;MFCCdn;MFCCdd2;eframesdn;eframesdd2];
MFCCfeats2(isnan(MFCCfeats2))=0;
MFCCfeats2_vector = MFCCfeats2(:);

E_distance = sqrt(sum((MFCCfeats2_vector - MFCCfeats1_vector).^2));

if E_distance <= s
    s = E_distance;
    n = i;
    m = AudioFile;
end

end
