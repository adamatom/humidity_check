import sys
import math


def dew_point(t_celsius, rh_percent):
    """
    Calculate the dew point temperature using the Magnus-Tetens approximation.

    Formula:
        dew_point = (b * alpha) / (a - alpha)
    where:
        alpha = ((a * t_celsius) / (b + t_celsius)) + ln(rh_percent / 100)

    Constants (typical):
        a = 17.27
        b = 237.7°C

    Parameters
    ----------
    t_celsius : float
        Indoor air temperature in °C.
    rh_percent : float
        Relative humidity in % (e.g., 50.0 for 50%).

    Returns
    -------
    float
        Dew point temperature in °C.
    """
    a = 17.27
    b = 237.7
    alpha = ((a * t_celsius) / (b + t_celsius)) + math.log(rh_percent / 100.0)
    dew_point_temp = (b * alpha) / (a - alpha)
    return dew_point_temp


def layer_temperatures(T_in, T_out, resistances):
    """
    Compute the temperature at each interface of a multi-layer assembly given the
    total indoor-outdoor temperature difference and the series thermal resistances.

    This is a 1D steady-state model. Each layer is characterized by an R-value (m²K/W).
    The temperature drop across each layer is proportional to its R-value relative
    to the total R.

    Parameters
    ----------
    T_in : float
        Indoor temperature (°C).
    T_out : float
        Outdoor temperature (°C).
    resistances : list of float
        A list of thermal resistances (R-values in m²K/W) in order from indoor to outdoor side.

    Returns
    -------
    list of float
        Temperatures at each interface, including T_in at the start and T_out at the end.
        For N layers, we get N+1 interface temperatures.
    """
    R_total = sum(resistances)
    delta_T = T_in - T_out
    temps = [T_in]
    current_T = T_in

    for R in resistances:
        # Temperature drop proportional to R of this layer
        dT = (R / R_total) * delta_T
        current_T = current_T - dT
        temps.append(current_T)
    return temps


def condensation_risk(T_in, rh, T_out, close_shutters=False):
    """
    Determine if there is a risk of condensation on the inner window surface.

    We consider a window assembly with layers:
      1. Indoor film resistance (R_in)
      2. Inner window pane (R_window1)
      3. Air gap (R_gap)
      4. Outer window pane (R_window2)
      5. Wooden shutters (R_shutter) if closed
      6. Outdoor film resistance (R_out)

    After computing the dew point from indoor temperature and relative humidity,
    we calculate the temperature profile across these layers. The critical point
    is the inner surface of the inner pane. If that temperature < dew point,
    condensation risk is True.

    Parameters
    ----------
    T_in : float
        Indoor temperature (°C).
    rh : float
        Indoor relative humidity (%).
    T_out : float
        Outdoor temperature (°C).
    close_shutters : bool, optional
        If True, shutters are closed, adding their resistance.

    Returns
    -------
    risk : bool
        True if condensation risk is present, otherwise False.
    dp : float
        Dew point temperature (°C).
    interface_temps : list of float
        Temperatures at each interface layer boundary.
    """
    dp = dew_point(T_in, rh)

    # Example R-values (m²K/W):
    # Film coefficients and pane approximations are from previous examples.
    R_in = 0.125      # Indoor film resistance
    R_window1 = 0.003 # Inner window pane (thin glass)
    R_gap = 0.5       # 10 cm gap with some convection (approx)
    R_window2 = 0.003 # Outer window pane
    R_out = 0.04      # Outdoor film resistance
    R_shutter = 0.26 if close_shutters else 0.0

    # Assemble the layer resistances
    resistances = [R_in, R_window1, R_gap, R_window2]
    if R_shutter > 0:
        resistances.append(R_shutter)
    resistances.append(R_out)

    interface_temps = layer_temperatures(T_in, T_out, resistances)
    inner_surface_temp = interface_temps[1]  # Just after indoor film, on inner pane surface

    risk = (inner_surface_temp < dp)
    return risk, dp, interface_temps, resistances


if __name__ == "__main__":
    # Command line arguments:
    # sys.argv[1] = Indoor Temperature (°C)
    # sys.argv[2] = Indoor Relative Humidity (%)
    # sys.argv[3] = Window Area (m²)
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <indoor_temperature_C> <indoor_relative_humidity_%> <window_area_m2>")
        sys.exit(1)

    try:
        T_in = float(sys.argv[1])
        rh = float(sys.argv[2])
        area = float(sys.argv[3])
    except ValueError:
        print("Error: All arguments must be numeric.")
        sys.exit(1)

    # Define a range of outside temperatures
    outside_temps = range(-10, 11, 5)  # -10, -5, 0, 5, 10 °C

    dp_current = dew_point(T_in, rh)
    print(f"\nIndoor Conditions: T_in={T_in:.1f}°C, RH={rh:.0f}%, Window Area={area:.2f} m²")
    print(f"Dew Point: {dp_current:.2f}°C\n")

    # Print header
    print("Outdoor Temp | Shutters Open (Risk, Heat Flow [W]) | Shutters Closed (Risk, Heat Flow [W])")
    for T_out in outside_temps:
        risk_open, dp_open, temps_open, R_open = condensation_risk(T_in, rh, T_out, close_shutters=False)
        risk_closed, dp_closed, temps_closed, R_closed = condensation_risk(T_in, rh, T_out, close_shutters=True)

        # Total R-values per unit area:
        R_total_open = sum(R_open)   # m²K/W
        R_total_closed = sum(R_closed)

        # Heat flow calculation (Q): Q = (T_in - T_out)*Area/R_total
        # This gives the steady-state heat flow in Watts.
        Q_open = (T_in - T_out)*area/R_total_open
        Q_closed = (T_in - T_out)*area/R_total_closed

        print(f"{T_out:>12}°C | {('Risk' if risk_open else 'No Risk'):>10}, {Q_open:8.2f} W "
              f"| {('Risk' if risk_closed else 'No Risk'):>13}, {Q_closed:8.2f} W")

