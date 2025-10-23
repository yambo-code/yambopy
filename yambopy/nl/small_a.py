import numpy as np
import math
from typing import Tuple, Optional, List, Union
import cmath

# External fields:
#
# SIN:      E(t)=\theta(t) * sin(t)
# SOFTSIN:  E(t)= (c*t^2 + b*t^3 + a*t^4 )* sin(t)  and 0 for t<=0
# DELTA:    E(t)= \delta(t)
# GAUSS:    E(t)= exp((t-t_o)^2/(2*sigma^2))
# THETA:    E(t)= \theta(t)
# RES:      E(t)= \theta(t)*exp(-i\omega t)
# ANTIRES:  E(t)= \theta(t)*exp(i\omega t)
# PULSE:    E(t)=(T-T_0 + sigma)**2 *(T-T_0 -sigma)**2/sigma**4 * cos(w t)
# QSSIN:    E(t)= exp((t-t_o)^2/(2*sigma^2))*sin(w t)
# SPULSE:   E(t)=((T-T_0)**2 - sigma**2)**2/sigma**4*cos(w t)
# PHHG:     E(t)=sin^2(\pi*(T-T_0)/sigma)*cos( w * t) 
# QSFIELD:  see below
# FROMFILE: shape of the electric field from file
#
# Linear frequency chirp coded for QSSIN pulse: 
# (see also https://www.rp-photonics.com/chirp.html)
#
#    s_sigma_chirp=sigma**2/(2.*(sigma**4+chirp**4))
#    c_sigma_chirp=chirp**2/(2.*(sigma**4+chirp**4))
#    sigma_eff=sqrt((sigma**4+chirp**4)/sigma**2)=sqrt(1./(2.*s_sigma_chirp))
#
#    sin(omega*(T-T0)+c_sigma_chirp*(T-T_0)**2) * exp(-(T-T_0)**2/(2.*sigma_eff**2) )
#
#   = [ exp( cI*(omega*(T-T0))*exp( cI*c_sigma_chirp*(T-T_0)**2) 
#      -exp(-cI*(omega*(T-T0))*exp(-cI*c_sigma_chirp*(T-T_0)**2) ]
#      * exp(-(T-T_0)**2 * s_sigma_chirp ) / 2.
#
#   = [ exp( cI*(omega*(T-T0))* exp( (s_sigma_chirp+cI*c_sigma_chirp) *(T-T_0)**2)
#     - exp(-cI*(omega*(T-T0))* exp( (s_sigma_chirp-cI*c_sigma_chirp) *(T-T_0)**2) ]
#

def small_a(T: float, dt: float, E_field: dict, order: int, envelop_only: Optional[bool] = False) -> Tuple[complex, complex]:
    """
    The vector potential is generally written as
    
    order=0  A (t)=-cEo  a (t) theta(t)
    order=1  A'(t)=-cEo (a'(t) theta(t)+a (t) delta(t))
    order=2  A"(t)=-cEo (a"(t) theta(t)+a'(t) delta(t)-a(t) sign(t))
    
    the functions theta,delta and sign can be the standard distributions
    or more fancy functions that can mimic the distributions.
    
    Note that A is evolved using A''(t) starting from A(0) and A'(0).
    """
    
    # Constants
    cI = 1j
    cONE = 1.0 + 0j
    cZERO = 0.0 + 0j
    pi = math.pi
    
    # Zeroing
    small_a_result = (cZERO, cZERO)
    f_t = [cZERO, cZERO]
    damp_func = 1.0
    
    envelop_only_ = envelop_only if envelop_only is not None else False
    
    # Field polarization
    if E_field.get("ef_pol") == "linear":
        n_fields = 1
    elif E_field.get("ef_pol") == "circular":
        n_fields = 2
    else:
        n_fields = 1  # default
    
    # Field parameters
    sigma = E_field.get("width", 0.0)
    chirp = E_field.get("chirp", 0.0)
    
    s_sigma_chirp = 0.0
    c_sigma_chirp = 0.0
    
    if abs(sigma) > 0.0 or abs(chirp) > 0.0:
        s_sigma_chirp = sigma**2 / (2.0 * (sigma**4 + chirp**4))
        c_sigma_chirp = chirp**2 / (2.0 * (sigma**4 + chirp**4))
    
    sigma_eff = math.sqrt((sigma**4 + chirp**4) / sigma**2) if sigma != 0 else 0.0
    
    # Parse field definitions
    field_defs = E_field.get("ef_name", "").split()
    
    Tloc = T
    if "RECT" in field_defs[0] and order == 0 and abs(T) >= sigma:
        Tloc = sigma
    
    # Determine field type and set parameters
    field_type = field_defs[0] if field_defs else ""
    
    W_0 = 0.0
    T_0 = 0.0
    damp_func = 1.0
    a, b, c = 0.0, 0.0, 0.0
    
    if field_type in ['STATIC', 'RECT', 'RECTSIN', 'SIN', 'DELTA']:
        # Fields which do not need T_0
        W_0 = 0.0
        T_0 = 0.0
        damp_func = 1.0
        
    elif field_type == 'FROM_FILE':
        # This would need file reading implementation
        i_file = get_field_file_index(field_defs[1])
        T_0 = field_from_file[0][0] * FS2AUT  # Assuming FS2AUT conversion
        W_0 = 0.0
        damp_func = 1.0
        
    elif field_type in ['SOFTSIN', 'THETA']:
        # Fields which do not need T_0 and with damp_func
        W_0 = 0.0
        T_0 = 0.0
        a = 3.0 / sigma**4 
        b = -8.0 / sigma**3
        c = 6.0 / sigma**2
        damp_func = 1.0
        if T < sigma and sigma > 0.0:
            damp_func = a * T**4 + b * T**3 + c * T**2
            
    elif field_type in ['GAUSS', 'QSSIN', 'QSFIELD', 'PULSE', 'SPULSE']:
        # Fields which need T_0
        W_0 = E_field.get("frequency", 0.0)
        T_0_fac = 3.0 * sigma_eff
        
        # Handle sigma specifications
        if "1SIGMA" in field_defs[1:] or "1SIGMA" in field_defs[2:]:
            T_0_fac = 1.0 * sigma_eff
        elif "2SIGMA" in field_defs[1:] or "2SIGMA" in field_defs[2:]:
            T_0_fac = 2.0 * sigma_eff
        elif "3SIGMA" in field_defs[1:] or "3SIGMA" in field_defs[2:]:
            T_0_fac = 3.0 * sigma_eff
        elif "4SIGMA" in field_defs[1:] or "4SIGMA" in field_defs[2:]:
            T_0_fac = 4.0 * sigma_eff
        elif "5SIGMA" in field_defs[1:] or "5SIGMA" in field_defs[2:]:
            T_0_fac = 5.0 * sigma_eff
            
        T_0 = pi / W_0 * (round(W_0 / pi * T_0_fac)) if W_0 != 0 else T_0_fac
        if "PULSE" in field_type:
            T_0 = T_0_fac
    
    # Initial and relative phases control
    fr_shift = [0.0, pi / 2.0]
    for i1 in range(len(field_defs)):
        if "PHPI180" in field_defs[i1]:
            fr_shift = [x + pi for x in fr_shift]  # 180 deg
        elif "PHPI120" in field_defs[i1]:
            fr_shift = [x + pi * 2.0 / 3.0 for x in fr_shift]  # 120 deg
        elif "PHPI90" in field_defs[i1]:
            fr_shift = [x + pi / 2.0 for x in fr_shift]  # 90 deg
        elif "PHPI60" in field_defs[i1]:
            fr_shift = [x + pi / 3.0 for x in fr_shift]  # 60 deg
        elif "PHPI30" in field_defs[i1]:
            fr_shift = [x + pi / 6.0 for x in fr_shift]  # 30 deg
        elif "PHPI20" in field_defs[i1]:
            fr_shift = [x + pi / 9.0 for x in fr_shift]  # 20 deg
    
    E_field["To"] = T_0
    
    for i_field in range(n_fields):
        W_field = E_field.get("frequency", 0.0)
        W_field_m1 = 1.0 / W_field if W_field > 0.0 else 0.0
        der_fac = W_field + 2.0 * c_sigma_chirp * (Tloc - T_0)
        
        # The frequency shift is applied in two cases
        # (i n_fields=2) to have a circular polarized pulse, and
        WtimesT = W_field * (Tloc - T_0) + fr_shift[i_field]
        # each frequency has a different initial phase
        if chirp > 0.0:
            WtimesT += c_sigma_chirp * (Tloc - T_0)**2
        
        if envelop_only_:
            f0t = cONE
            f1t = cONE
        else:
            # CONTROL RES / ANTIRES case
            cos_wt = math.cos(WtimesT)
            sin_wt = math.sin(WtimesT)
            exp_iwt = complex(cos_wt, sin_wt)
            
            f0t = complex(cos_wt, 0.0)
            f1t = complex(sin_wt, 0.0)
            
            if len(field_defs) > 1 and field_defs[1] == "ANTIRES":
                f0t = 0.5 * exp_iwt
                f1t = -cI * 0.5 * exp_iwt
            elif len(field_defs) > 1 and field_defs[1] == "RES":
                f0t = 0.5 * exp_iwt.conjugate()
                f1t = cI * 0.5 * exp_iwt.conjugate()
        
        EXPf = math.exp(-(T - T_0)**2 / (2.0 * sigma_eff**2)) if sigma_eff != 0 else 0.0
        
        # Field type specific calculations
        f_now = cZERO
        
        if field_type == 'FROM_FILE':
            # File reading implementation needed
            i_T = int(round((T - T_0) / (dt / 2.0))) + 1
            if i_T <= 0 or envelop_only_:
                f_now = 0
            else:
                # This would need actual file data access
                pass
                
        elif field_type == 'STATIC':
            if order == 0:
                f_now = T
            elif order == 1:
                f_now = 1.0
            elif order == 2:
                f_now = 0.0
                
        elif field_type == 'RECT':
            if order == 0:
                f_now = Tloc
            elif order == 1:
                f_now = theta_function(sigma - T, dt, 0)  # theta function
            elif order == 2:
                f_now = -theta_function(sigma - T, dt, 1)  # delta function
                
        elif field_type == 'RECTSIN':
            if chirp > 0.0:
                raise ValueError("chirp not implemented with " + field_type)
                
            if order == 0:
                f_now = -(f0t - 1.0) * W_field_m1
            elif order == 1:
                f_now = theta_function(sigma - T, dt, 0) * f1t  # theta function
            elif order == 2:
                f_now = (-theta_function(sigma - T, dt, 1) * f1t +  # delta function
                         theta_function(sigma - T, dt, 0) * f0t * der_fac)
                        
        elif field_type == 'SIN':
            if chirp > 0.0:
                raise ValueError("chirp not implemented with " + field_type)
                
            if order == 0:
                f_now = -damp_func * (f0t - 1.0) * W_field_m1
            elif order == 1:
                f_now = damp_func * f1t
            elif order == 2:
                f_now = damp_func * f0t * der_fac
                
        elif field_type == 'SOFTSIN':
            if chirp > 0.0:
                raise ValueError("chirp not implemented with " + field_type)
                
            if order == -1:
                f_now = -2.0
            elif order == 0:
                f_now = -damp_func * (f0t - 1.0) * W_field_m1
            elif order == 1:
                f_now = damp_func * f1t
            elif order == 2:
                f_now = damp_func * f0t * der_fac
                
        elif field_type == 'THETA':
            if order == 0:
                f_now = damp_func * T
            elif order == 1:
                f_now = damp_func
            elif order == 2:
                f_now = 0.0
                
        elif field_type == 'DELTA':
            if order == -1:
                f_now = 1.0
            elif order == 0:
                f_now = 1.0
            elif order > 0:
                f_now = 0.0
                
        elif field_type == 'PHHG':
            if chirp > 0.0:
                raise ValueError("chirp not implemented with " + field_type)
                
            sarg = pi * (T - T_0) / sigma
            WT = W_field * T
            
            if (T - T_0 <= 0.0) or (T - T_0 >= sigma and order > 0):
                f_now = 0.0
            elif T - T_0 >= sigma and order == 0:
                Tl = sigma + T_0
                WT = W_field * Tl
                f_now = (-(sigma * math.sin(((sigma * W_field + 2 * pi) * Tl - 2 * pi * T_0) / sigma)) / 
                        (4 * (sigma * W_field + 2 * pi)) -
                        (sigma * math.sin(((sigma * W_field - 2 * pi) * Tl + 2 * pi * T_0) / sigma)) / 
                        (4 * (sigma * W_field - 2 * pi)) + math.sin(WT) / (2 * W_field))
            else:
                if order == 0:
                    f_now = (-(sigma * math.sin(((sigma * W_field + 2 * pi) * T - 2 * pi * T_0) / sigma)) / 
                            (4 * (sigma * W_field + 2 * pi)) -
                            (sigma * math.sin(((sigma * W_field - 2 * pi) * T + 2 * pi * T_0) / sigma)) / 
                            (4 * (sigma * W_field - 2 * pi)) + math.sin(WT) / (2 * W_field))
                elif order == 1:
                    f_now = math.sin(sarg)**2 * math.cos(WT)
                elif order == 2:
                    f_now = ((2 * pi * math.cos(WT) * math.cos(sarg) * math.sin(sarg)) / sigma - 
                            W_field * math.sin(WT) * math.sin(sarg)**2)
                            
        elif field_type == 'GAUSS':
            if order == 0:
                f_now = (sigma_eff * math.sqrt(pi / 2.0) * 
                        (math.erf((T - T_0) / (sigma_eff * math.sqrt(2.0))) + 1.0))
            elif order == 1:
                f_now = EXPf
            elif order == 2:
                f_now = -EXPf * (T - T_0) / sigma_eff**2
                
        elif field_type == 'QSSIN':
            # This would need FADEVA function implementation
            # Placeholder implementation
            if order == 0:
                f_now = cZERO  # Placeholder
            elif order == 1:
                f_now = f1t * EXPf
            elif order == 2:
                f_now = (der_fac * f0t - (T - T_0) * f1t / sigma_eff**2) * EXPf
                
        elif field_type == 'QSFIELD':
            if order == 0:
                f_now = f1t * EXPf
            elif order == 1:
                f_now = (der_fac * f0t - (T - T_0) * f1t / sigma**2) * EXPf
            elif order == 2:
                f_now = ((-der_fac * f1t - f1t / sigma**2 - 
                         der_fac * (T - T_0) * f1t / sigma**2 - 
                         (T - T_0) * (der_fac * f0t - (T - T_0) * f1t / sigma**2) / sigma**2) * EXPf)
                f_now = f_now / der_fac
                
        elif field_type == 'PULSE':
            if chirp > 0.0:
                raise ValueError("chirp not implemented with " + field_type)
                
            if abs(T - T_0) < sigma:
                if order == 0:
                    f_now = 0.0
                elif order == 1:
                    f_now = ((T - T_0 + sigma)**2 * (T - T_0 - sigma)**2 / sigma**4 * f0t)
                elif order == 2:
                    f_now = (4.0 * (T - T_0 + sigma) * (T - T_0 - sigma)**2 / sigma**4 * f0t -
                            (T - T_0 + sigma)**2 * (T - T_0 - sigma)**2 / sigma**4 * W_field * f1t)
                            
        elif field_type == 'SPULSE':
            if chirp > 0.0:
                raise ValueError("chirp not implemented with " + field_type)
                
            T_0 = sigma
            W_0 = W_field
            f_now = cZERO
            
            if abs(T - T_0) < sigma:
                # Simplified implementation - full Fortran translation would be more complex
                if order == 1:
                    f_now = ((T - T_0)**2 - sigma**2)**2 / sigma**4 * f0t
                    
        f_t[i_field] = f_t[i_field] + f_now
    
    return tuple(f_t)

# Helper functions that would need implementation
def theta_function(x: float, dt: float, order: int) -> float:
    """
    Placeholder for theta function implementation
    """
    # This would need proper implementation based on the Fortran version
    if order == 0:
        return 1.0 if x > 0 else 0.0
    elif order == 1:
        # Approximate delta function
        return 1.0 / dt if abs(x) < dt/2 else 0.0
    return 0.0

def get_field_file_index(filename: str) -> int:
    """
    Placeholder for file index retrieval
    """
    return 0

# Global variables that would need to be defined
FS2AUT = 1.0  # Placeholder - actual conversion factor needed
field_from_file = []  # Placeholder for file data
