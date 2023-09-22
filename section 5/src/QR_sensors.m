function [P,pivots] = QR_sensors(Vr,p)

[n,r] = size(Vr);
p = min([p,n]); p = max([p,r]);     % r <= p <= n
I = speye(n);
if p==r                             % QR sensor selection, p=r (Q-DEIM)
[~,~,pivots] = qr(Vr','vector');
elseif p>r                          % Oversampled QR sensors, p>r
[~,~,pivots] = qr(Vr*Vr','vector');
end
pivots = pivots(1:p);
P = I(:,pivots);
end