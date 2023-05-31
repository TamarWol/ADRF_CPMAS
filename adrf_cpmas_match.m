% PROSPR sequence

function mag=adrf_cpmas_match(spin_system,parameters,H,R,K)

% Compose Liouvillian
L=H+1i*R+1i*K;

% Get pulse operators 
Hp=operator(spin_system,'L+',parameters.spins{1});
Hx=(Hp+Hp')/2;  Hx=kron(speye(parameters.spc_dim),Hx);
Hy=(Hp-Hp')/2i; Hy=kron(speye(parameters.spc_dim),Hy); %#ok<NASGU>
Cp=operator(spin_system,'L+',parameters.spins{2});
Cx=(Cp+Cp')/2;  Cx=kron(speye(parameters.spc_dim),Cx);
Cy=(Cp-Cp')/2i; Cy=kron(speye(parameters.spc_dim),Cy); %#ok<NASGU>

% Apply the ADRF
rho_adrf=shaped_pulse_xy(spin_system,L,{Hx},{fliplr(parameters.ramp_amp)},...
                         parameters.ramp_dt*ones(size(parameters.ramp_amp)),...
                         parameters.rho0);

% Preallocate answer
mag=zeros(size(parameters.cw_amp));

% Apply the spin-lock period
for n=1:numel(parameters.cw_amp)

    % Make spin-lock operator
    L_SL=L+parameters.cw_amp(n)*Cx;
    
    % Run the evolution
    rho=evolution(spin_system,L_SL,[],rho_adrf,parameters.cw_dur,1,'final');

    % Get the Cx observable
    mag(n)=real(trace(parameters.coil'*rho));
    
end

end

