% cy figure2 25/9/5

tic;
clear;
set(groot, 'DefaultAxesFontSize', 20); 
set(groot, 'DefaultAxesFontName', 'Times New Roman');  
set(groot, 'DefaultTextFontName', 'Times New Roman'); 

distance = 3;   % normalized to dispersion length

delta3= 0.0;   % beta3/(6*T0*abs(beta2))
delta4 = 0.0;   % beta4/(24*T0^2*abs(beta2))
chirp0 = 0;     % input pulse chirp (default value)
s=0;
T0 = 50;            % pulse width in fs
sbeta2 = -1;        % sign of beta_2

nt = 2^15;          % FFT points
Tmax = 4;          % FFT window size
step_num = 5000;       % We only need one step forward and one step backward
zstep = step_num;
dtau = (2 * Tmax) / nt;        % step size in tau
%--- Temporal grid and frequency array ---
tau = (-nt/2:nt/2-1) * dtau;  % temporal grid
omega = fftshift((pi / Tmax) * (-nt/2:nt/2-1));  % frequency array

deltaz = distance / step_num;  % step size in z
ng = step_num/zstep;        % store data for 3D.

% mps=0:1:10;
mps = 0;
m_p=0;
E_collect=nan(length(mps),step_num+1);
i_m_ps=1;

for mshape = 2
    for N=3:4
        %%
        
        uu_0_ref = exp(-0.5*(1+1i*chirp0).*tau.^(2*mshape));  % super-Gaussian
        uu = uu_0_ref;
        uu_0_no_noise = zeros(step_num, nt);
        uu_0_no_noise(1,:) = uu;

        fR = 0;fb=0;
        tau1=12.2/T0; tau2=32/T0; tau_b=96/T0;
        h = (1-fb)*(tau1^2 + tau2^2)/(tau1*tau2^2)*exp(-tau/tau2).*sin(tau/tau1)...
            +fb*((2*tau_b-tau)/tau_b^2).*exp(-tau/tau_b);
        h(1:nt/2) = 0;      % causality

        dispersion = exp(1i*deltaz*(0.5*sbeta2*omega.^2 +...
            delta3*omega.^3 + delta4*omega.^4));
        hhz = 1i*N^2*deltaz;


        temp = uu .* exp(abs(uu).^2 .* hhz / 2);
        for n = 1:step_num
            f_temp = ifft(temp) .* dispersion;
            uu = fft(f_temp);

            sst1 = deriv1(abs(uu).^2,dtau);
            sst2 = deriv1(uu,dtau);
            sst = 1i*(1-fR)*s*(sst1 + conj(uu).*sst2);
            P = uu.*conj(uu);
            convl = (nt*dtau)*fft(ifft(h).*ifft(P));
            convl = fftshift(convl);
            sst3 = deriv1(convl,dtau);
            sstnew = 1i*s*fR*(sst3 + (conj(uu).*convl)./(P+eps).*sst2);

            temp = uu.*exp(((1-fR)*(abs(uu).^2)+ sst+ fR*convl+ sstnew).*hhz);
            % temp = temp .* exp(1i * potential * deltaz);

            if mod(n,ng) == 0
                uu_0_no_noise((n/ng+1), :) = temp;
            end
        end

        % 能量计算
        U_ref = uu_0_no_noise;p=3;lamda=-N^2;
        dU = nan(size(U_ref));
        for i = 1:size(U_ref,1)
            dU(i,:) = gradient(U_ref(i,:));
        end
        E_xi = .5*sum(abs(dU).^2,2)*dtau + ...
            2*(-N^2)/(p+1)*sum(abs(U_ref).^(p+1),2)*dtau;
        term3 = sum((0.5*m_p*tau.^2) .* abs(U_ref).^2, 2) * dtau;
        signal0 = E_xi + term3;  
        % ======  ======

        sigma = 0.1;
        uu_0 = exp(-0.5*(1+1i*chirp0).*tau.^(2*mshape)); %super-Gaussian
        % NOISE
        noise = sigma * randn(size(tau)); 
        uu_noisy = uu_0 + noise; 
        uu_0=uu_noisy;
        uu_initial = uu_0;
        %% 
        % Definition of Raman response function
        fR = 0;fb=0;
        tau1=12.2/T0; tau2=32/T0; tau_b=96/T0;
        h = (1-fb)*(tau1^2 + tau2^2)/(tau1*tau2^2)*exp(-tau/tau2).*sin(tau/tau1)...
            +fb*((2*tau_b-tau)/tau_b^2).*exp(-tau/tau_b);
        h(1:nt/2) = 0;      % causality
        %---Input Field profile
        uu = uu_0;
        % uu = uu_noisy;
        %-----Input 3D data
        uu_0=zeros(step_num,nt);
        uu_0(1, :)=uu;
        %----store dispersive phase shifts to speedup code
        dispersion = exp(1i*deltaz*(0.5*sbeta2*omega.^2 +...
            delta3*omega.^3 + delta4*omega.^4));
        hhz = 1i*N^2*deltaz;

        % 
        potential = -0.5 * m_p * tau.^2; 

        %*********[ Beginning of MAIN Loop]***********
        temp = uu.*exp(abs(uu).^2.*hhz/2);    % 
        for n = 1:step_num
            f_temp = ifft(temp).*dispersion; 
            uu = fft(f_temp);

            sst1 = deriv1(abs(uu).^2,dtau);
            sst2 = deriv1(uu,dtau);
            sst = 1i*(1-fR)*s*(sst1 + conj(uu).*sst2);
            P = uu.*conj(uu);
            convl = (nt*dtau)*fft(ifft(h).*ifft(P));
            convl = fftshift(convl);
            sst3 = deriv1(convl,dtau);
            sstnew = 1i*s*fR*(sst3 + (conj(uu).*convl)./(P+eps).*sst2);

            temp = uu.*exp(((1-fR)*(abs(uu).^2)+ sst+ fR*convl+ sstnew).*hhz);

            temp = temp .* exp(1i * potential * deltaz);

            spect = fftshift(ifft(temp));
            if mod(n,ng) == 0
                uu_0((n/ng+1),:) = temp;
            end
        end
        %% 
        z=linspace(0,distance,zstep+1); % normalised distance

        U=uu_0;p=3;lamda=-N^2;
        
        dU = nan(size(U));
        
        for i = 1:size(U,1)
            
            dU(i, :) = gradient(U(i, :)); 
        end

        E_xi = .5*sum(abs(dU).^2,2)*dtau + ...
            2*lamda/(p+1)*sum(abs(U).^(p+1),2)*dtau;

        % ∫ τ^2/2 * ψ^2 dτ
        term3 = sum((0.5*m_p*tau.^2) .* abs(U).^2, 2) * dtau;
        E_collect(i_m_ps,:)=E_xi+term3;
        signal=E_xi+term3;
        %% 
        figure
        subplot(1,4,1)
        plot(tau,uu_initial)
        xlabel('Distance (z/L_D)');ylabel('Intensity')

        % title(sprintf('Initial pluse'))
        subplot(1,4,2)
        imagesc(z,tau,abs(uu_0).^2');hold on
        set(gca,'XDir','normal')
        % ylim([-x_uu_length x_uu_length]);
        xlabel('Distance (z/L_D)'); ylabel('Time (t/T_0)');
        colorbar
        set(gcf,Position=[30.6000 242 1.3736e+03 450.2000])
        % title('Intensity-z')

        subplot(1,4,3)
        plot(z, signal, 'b', 'DisplayName', 'With noise');
        hold on
        plot(z, signal0, 'r--', 'DisplayName', 'No noise');
        xlabel('Distance (z/L_D)');
        ylabel('Energy');
        % title('Energy vs z');
        grid on;
        % legend;

        signal = signal - mean(signal); 
        
        padding_factor = 1;
        window = hamming(length(signal))';
        windowed_signal = signal.* window;

        
        % N = length(signal);
        N_fft = padding_factor * step_num;

                fft_signal = abs(fftshift(fft(windowed_signal, N_fft)));
        freq_axis = (-N_fft/2:N_fft/2-1) * (1 / (distance * padding_factor));

        
        positive_indices = freq_axis >= 0;
        fft_signal_pos = fft_signal(positive_indices);
        freq_axis_pos = freq_axis(positive_indices);

        
        [pks_peaks, locs_peaks] = findpeaks(fft_signal_pos, freq_axis_pos);
        peak_median = median(pks_peaks);

        
        valid_indices = (locs_peaks >= 0 & locs_peaks <= 6);
        valid_pks_peaks = pks_peaks(valid_indices);
        valid_locs_peaks = locs_peaks(valid_indices);
        % 
        % figure
        subplot(1,4,4)
        plot(freq_axis_pos, fft_signal_pos);
        hold on
        plot(valid_locs_peaks, valid_pks_peaks, 'ro') 
        hold off
        xlim([0 15])
        xlabel('Wavevector k')
        ylabel('FFT Magnitude')
        % title('FFT Spectrum')
        grid on
       
        title_str = sprintf("cyfig2N=%.1f;m=%.1f;noise=%.2f",N,mshape,sigma);
        % sgtitle(title_str);

        set(gcf,"Position", [195.4000 284.2000 1312 264.8000])
        
        saveas(gcf,strcat(title_str,'.png'));
        saveas(gcf,strcat(title_str,'.fig'));

        % close(gcf)

      

    end
end

%%

toc
function df1=deriv1(y,h)
% first-order derivative of y array accurate to h^2
% derivatives at end points are from page 527 of Chapra's book
% assumes equispaced x array with step size h
n=length(y);
df1 = (y(3:end) - y(1:end-2)); 	% central differences
d1= -y(3)+4*y(2)-3*y(1);		% first point
dn= 3*y(n)-4*y(n-1)+y(n-2);		% last point
df1 = [d1 df1 dn]./(2*h);		% derivative at all points
end