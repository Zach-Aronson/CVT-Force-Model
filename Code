
import csv
import json
import math
import webbrowser
from pathlib import Path

import plotly.graph_objects as go


# =============================================================================
# DEFAULT SETTINGS — EDIT THESE TO CHANGE THE PAGE'S STARTING VALUES
# =============================================================================

# --- Grid and graph settings -------------------------------------------------

MIN_RPM = 0
MAX_RPM = 5000
RPM_STEP = 100

POSITION_POINTS = 25

TORQUE_FEEDBACK_MIN_LBF = 0.0
TORQUE_FEEDBACK_MAX_LBF = 200.0
TORQUE_FEEDBACK_STEP_LBF = 5.0

POINT_SIZE = 2.8
POINT_OPACITY = 0.70
USE_CUBE_ASPECT_RATIO = True

# Fixed color range. Values outside the range remain in the calculation
# and hover text, but their displayed color is saturated.
NET_FORCE_COLOR_MIN_LBF = -100.0
NET_FORCE_COLOR_MAX_LBF = 100.0

DEFAULT_VIEW_MODE = "full"       # "full", "slices", or "cutaway"
DEFAULT_RPM_CUT = 3000
DEFAULT_SHIFT_CUT_IN = 0.438
DEFAULT_TORQUE_FEEDBACK_CUT = 100.0

OUTPUT_FILE_NAME = "cvt_torque_feedback_cube.html"
DEFAULT_CSV_FILE_NAME = "cvt_torque_feedback_cube.csv"
SAVE_DEFAULT_CSV = True


# --- Primary clutch defaults -------------------------------------------------

PRIMARY_SPRING_RATE = 60.0       # lbf/in
PRIMARY_PRETENSION = 2.20        # in
PRIMARY_WEIGHTS = 1.144          # lbf total, all flyweights combined

# Renamed from Primary Weight Leverage. This remains an effective,
# calibrated radius in the current Excel-style force model.
PRIMARY_WEIGHT_RADIUS_IN = 2.75  # in

PRIMARY_RAMP_ANGLE_DEG = 28.55   # deg
PRIMARY_RAMP_REFERENCE_DEG = 28.55

# This is the calibrated force factor from the original spreadsheet:
# -SIN(28.55) / 12
PRIMARY_FORCE_FACTOR = 0.02268233


# --- Secondary clutch defaults -----------------------------------------------

SECONDARY_SPRING_RATE = 20.0     # lbf/in
SECONDARY_PRETENSION = 2.00      # in

# Stored for future physical torque-cam work. In this parameter sweep,
# torque feedback is the independent Z-axis, so these do not drive the cube.
SECONDARY_CAM_ANGLE_DEG = 32.0
SECONDARY_CAM_DIAMETER_IN = 3.50
SECONDARY_CAM_RADIUS_IN = SECONDARY_CAM_DIAMETER_IN / 2.0
SECONDARY_CAM_CALIBRATION = 0.20


# --- CVT / belt geometry -----------------------------------------------------

BELT_WIDTH = 0.875               # in
LOW_RATIO = 3.90
HIGH_RATIO = 0.90


# --- CH440 dyno data ---------------------------------------------------------
# Format: (RPM, torque in ft-lb)

ENGINE_TORQUE_CURVE = [
    (2400, 18.5),
    (2600, 18.1),
    (2800, 17.4),
    (3000, 16.6),
    (3200, 15.4),
    (3400, 14.5),
    (3600, 13.5),
]


# =============================================================================
# ENGINE TORQUE CURVE
# =============================================================================

def solve_linear_system(matrix, vector):
    """Solve a square linear system using Gaussian elimination."""

    size = len(vector)
    augmented = [
        [float(matrix[row][column]) for column in range(size)]
        + [float(vector[row])]
        for row in range(size)
    ]

    for pivot_column in range(size):
        pivot_row = max(
            range(pivot_column, size),
            key=lambda row: abs(augmented[row][pivot_column]),
        )

        if abs(augmented[pivot_row][pivot_column]) < 1e-18:
            raise ValueError("Torque-fit matrix is singular.")

        augmented[pivot_column], augmented[pivot_row] = (
            augmented[pivot_row],
            augmented[pivot_column],
        )

        pivot_value = augmented[pivot_column][pivot_column]

        for column in range(pivot_column, size + 1):
            augmented[pivot_column][column] /= pivot_value

        for row in range(size):
            if row == pivot_column:
                continue

            factor = augmented[row][pivot_column]

            for column in range(pivot_column, size + 1):
                augmented[row][column] -= (
                    factor * augmented[pivot_column][column]
                )

    return [augmented[row][size] for row in range(size)]


def fit_constrained_cubic_through_origin(data_points):
    """
    Fit a least-squares curve:

        torque = a*RPM^3 + b*RPM^2 + c*RPM

    It is constrained through (0 RPM, 0 ft-lb). It is a smooth best-fit
    curve, so it approximates rather than exactly passes through every
    measured dyno point. This gives a more reasonable extrapolation to
    5000 RPM than a constrained quadratic.
    """

    sum_x2 = 0.0
    sum_x3 = 0.0
    sum_x4 = 0.0
    sum_x5 = 0.0
    sum_x6 = 0.0

    sum_xy = 0.0
    sum_x2y = 0.0
    sum_x3y = 0.0

    for rpm, torque in data_points:
        sum_x2 += rpm**2
        sum_x3 += rpm**3
        sum_x4 += rpm**4
        sum_x5 += rpm**5
        sum_x6 += rpm**6

        sum_xy += rpm * torque
        sum_x2y += (rpm**2) * torque
        sum_x3y += (rpm**3) * torque

    matrix = [
        [sum_x6, sum_x5, sum_x4],
        [sum_x5, sum_x4, sum_x3],
        [sum_x4, sum_x3, sum_x2],
    ]

    vector = [
        sum_x3y,
        sum_x2y,
        sum_xy,
    ]

    coefficient_a, coefficient_b, coefficient_c = solve_linear_system(
        matrix,
        vector,
    )

    return coefficient_a, coefficient_b, coefficient_c


TORQUE_FIT_A, TORQUE_FIT_B, TORQUE_FIT_C = (
    fit_constrained_cubic_through_origin(ENGINE_TORQUE_CURVE)
)


def engine_torque_at_rpm(rpm):
    """
    Calculate engine torque from the constrained cubic best-fit curve.

    Torque is clamped to zero if extrapolation ever predicts a negative
    value. Values from 0–2400 RPM and 3600–5000 RPM are extrapolations,
    not measured dyno data.
    """

    torque = (
        TORQUE_FIT_A * rpm**3
        + TORQUE_FIT_B * rpm**2
        + TORQUE_FIT_C * rpm
    )

    return max(0.0, torque)


# =============================================================================
# CVT MODEL
# =============================================================================

def clamp(value, minimum, maximum):
    """Clamp a value inside a minimum and maximum range."""
    return max(minimum, min(value, maximum))


def nearest_index(values, target):
    """Return index of the grid point nearest to target."""
    return min(
        range(len(values)),
        key=lambda index: abs(values[index] - target),
    )


def cvt_ratio(shift_position):
    """Linearly map shift position from low ratio to high ratio."""

    if BELT_WIDTH <= 0:
        raise ValueError("BELT_WIDTH must be greater than zero.")

    fraction = clamp(shift_position / BELT_WIDTH, 0.0, 1.0)

    return LOW_RATIO + fraction * (HIGH_RATIO - LOW_RATIO)


def primary_ramp_multiplier(ramp_angle_deg):
    """
    Scale the original calibrated primary force factor by ramp angle.

    At 28.55 deg, the multiplier is exactly 1.0, preserving the current
    force model. This uses a sine relationship as a simple approximation.
    """

    reference_sine = math.sin(
        math.radians(PRIMARY_RAMP_REFERENCE_DEG)
    )

    if abs(reference_sine) < 1e-12:
        return 1.0

    return math.sin(math.radians(ramp_angle_deg)) / reference_sine


def primary_force(
    rpm,
    shift_position,
    primary_spring_rate,
    primary_pretension,
    primary_weights,
    primary_weight_radius_in,
    primary_ramp_angle_deg,
):
    """
    Calculate primary applied force.

    Flyweight force =
        RPM * total flyweight weight * effective weight radius
        * calibrated force factor * ramp multiplier

    Primary force =
        flyweight force - primary spring force
    """

    flyweight_force = (
        rpm
        * primary_weights
        * primary_weight_radius_in
        * PRIMARY_FORCE_FACTOR
        * primary_ramp_multiplier(primary_ramp_angle_deg)
    )

    primary_spring_force = (
        primary_spring_rate
        * (primary_pretension - shift_position)
    )

    return flyweight_force - primary_spring_force


def secondary_spring_force(
    shift_position,
    secondary_spring_rate,
    secondary_pretension,
):
    """Calculate secondary compression-spring force."""

    compression = secondary_pretension + shift_position

    return secondary_spring_rate * compression


def calculate_point(
    rpm,
    shift_position,
    torque_feedback_lbf,
    primary_spring_rate=PRIMARY_SPRING_RATE,
    primary_pretension=PRIMARY_PRETENSION,
    secondary_spring_rate=SECONDARY_SPRING_RATE,
    secondary_pretension=SECONDARY_PRETENSION,
    primary_weights=PRIMARY_WEIGHTS,
    primary_weight_radius_in=PRIMARY_WEIGHT_RADIUS_IN,
    primary_ramp_angle_deg=PRIMARY_RAMP_ANGLE_DEG,
):
    """
    Calculate one point in the 3D parameter sweep.

    Net shift force =
        primary force
        - secondary spring force
        - torque feedback
    """

    ratio = cvt_ratio(shift_position)
    engine_torque = engine_torque_at_rpm(rpm)

    primary = primary_force(
        rpm,
        shift_position,
        primary_spring_rate,
        primary_pretension,
        primary_weights,
        primary_weight_radius_in,
        primary_ramp_angle_deg,
    )

    secondary_spring = secondary_spring_force(
        shift_position,
        secondary_spring_rate,
        secondary_pretension,
    )

    secondary_total = secondary_spring + torque_feedback_lbf
    net_force = primary - secondary_total

    return {
        "rpm": rpm,
        "shift_position_in": shift_position,
        "torque_feedback_lbf": torque_feedback_lbf,
        "cvt_ratio": ratio,
        "engine_torque_ftlb": engine_torque,
        "primary_force_lbf": primary,
        "secondary_spring_lbf": secondary_spring,
        "secondary_total_lbf": secondary_total,
        "net_shift_force_lbf": net_force,
    }


# =============================================================================
# GRID / CSV HELPERS
# =============================================================================

def make_rpm_values():
    """Build RPM grid and retain the original measured dyno RPM values."""

    if RPM_STEP <= 0:
        raise ValueError("RPM_STEP must be greater than zero.")

    values = set(range(MIN_RPM, MAX_RPM + 1, RPM_STEP))
    values.add(MIN_RPM)
    values.add(MAX_RPM)

    for dyno_rpm, _ in ENGINE_TORQUE_CURVE:
        if MIN_RPM <= dyno_rpm <= MAX_RPM:
            values.add(dyno_rpm)

    return sorted(values)


def make_shift_positions():
    """Build equally spaced shift positions from 0 to BELT_WIDTH."""

    if POSITION_POINTS < 2:
        raise ValueError("POSITION_POINTS must be at least 2.")

    return [
        BELT_WIDTH * index / (POSITION_POINTS - 1)
        for index in range(POSITION_POINTS)
    ]


def make_torque_feedback_values():
    """Build the torque-feedback Z-axis grid."""

    if TORQUE_FEEDBACK_STEP_LBF <= 0:
        raise ValueError(
            "TORQUE_FEEDBACK_STEP_LBF must be greater than zero."
        )

    if TORQUE_FEEDBACK_MAX_LBF < TORQUE_FEEDBACK_MIN_LBF:
        raise ValueError(
            "TORQUE_FEEDBACK_MAX_LBF must be greater than or equal to "
            "TORQUE_FEEDBACK_MIN_LBF."
        )

    values = []
    current = TORQUE_FEEDBACK_MIN_LBF

    while current <= TORQUE_FEEDBACK_MAX_LBF + 1e-9:
        values.append(round(current, 6))
        current += TORQUE_FEEDBACK_STEP_LBF

    if values[-1] != TORQUE_FEEDBACK_MAX_LBF:
        values.append(TORQUE_FEEDBACK_MAX_LBF)

    return values


def write_csv(rows):
    """Write the default calculation point cloud to a CSV file."""

    output_path = Path(DEFAULT_CSV_FILE_NAME).resolve()

    fieldnames = [
        "rpm",
        "shift_position_in",
        "torque_feedback_lbf",
        "cvt_ratio",
        "engine_torque_ftlb",
        "primary_force_lbf",
        "secondary_spring_lbf",
        "secondary_total_lbf",
        "net_shift_force_lbf",
    ]

    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return output_path


# =============================================================================
# HTML CONTROL BUILDERS
# =============================================================================

def model_input(input_id, label, value, minimum, maximum, step, units):
    """Return one numeric control for the HTML sidebar."""

    return f"""
    <div class="model-row">
      <label for="{input_id}">{label}</label>
      <div class="input-group">
        <input id="{input_id}" type="number"
               value="{value}" min="{minimum}" max="{maximum}"
               step="{step}">
        <span>{units}</span>
      </div>
    </div>
    """


def make_graph():
    """Create the interactive HTML point cloud and sidebar controls."""

    rpm_values = make_rpm_values()
    shift_values = make_shift_positions()
    torque_values = make_torque_feedback_values()

    points_rpm = []
    points_shift = []
    points_torque = []
    default_colors = []
    default_customdata = []
    default_rows = []

    for torque_feedback in torque_values:
        for shift_position in shift_values:
            for rpm in rpm_values:
                point = calculate_point(
                    rpm,
                    shift_position,
                    torque_feedback,
                )

                points_rpm.append(point["rpm"])
                points_shift.append(point["shift_position_in"])
                points_torque.append(point["torque_feedback_lbf"])

                default_colors.append(point["net_shift_force_lbf"])

                default_customdata.append([
                    point["cvt_ratio"],
                    point["engine_torque_ftlb"],
                    point["primary_force_lbf"],
                    point["secondary_spring_lbf"],
                    point["secondary_total_lbf"],
                    point["net_shift_force_lbf"],
                ])

                default_rows.append(point)

    rpm_cut_index = nearest_index(rpm_values, DEFAULT_RPM_CUT)
    shift_cut_index = nearest_index(shift_values, DEFAULT_SHIFT_CUT_IN)
    torque_cut_index = nearest_index(
        torque_values,
        DEFAULT_TORQUE_FEEDBACK_CUT,
    )

    figure = go.Figure(
        data=[
            go.Scatter3d(
                x=points_rpm,
                y=points_shift,
                z=points_torque,
                mode="markers",
                marker=dict(
                    size=POINT_SIZE,
                    opacity=POINT_OPACITY,
                    color=default_colors,
                    colorscale="RdBu_r",
                    cmin=NET_FORCE_COLOR_MIN_LBF,
                    cmax=NET_FORCE_COLOR_MAX_LBF,
                    colorbar=dict(
                        title="Net Shift Force (lbf)",
                        tickvals=[-100, -50, 0, 50, 100],
                    ),
                ),
                customdata=default_customdata,
                hovertemplate=(
                    "RPM: %{x:.0f}<br>"
                    "Shift Position: %{y:.3f} in<br>"
                    "Torque Feedback: %{z:.1f} lbf<br>"
                    "CVT Ratio: %{customdata[0]:.3f}<br>"
                    "Engine Torque: %{customdata[1]:.2f} ft-lb<br>"
                    "Primary Force: %{customdata[2]:.1f} lbf<br>"
                    "Secondary Spring: %{customdata[3]:.1f} lbf<br>"
                    "Secondary Total: %{customdata[4]:.1f} lbf<br>"
                    "Net Shift Force: %{customdata[5]:.1f} lbf"
                    "<extra></extra>"
                ),
                showlegend=False,
            )
        ]
    )

    aspect_mode = "cube" if USE_CUBE_ASPECT_RATIO else "auto"

    figure.update_layout(
        title=(
            "CVT Torque Feedback Parameter Sweep"
            "<br><sup>Point Height = Torque Feedback | "
            "Point Color = Net Shift Force | "
            "Color Range Fixed at −100 to +100 lbf</sup>"
        ),
        scene=dict(
            xaxis=dict(title="Engine RPM", range=[MIN_RPM, MAX_RPM]),
            yaxis=dict(title="Shift Position (in)", range=[0, BELT_WIDTH]),
            zaxis=dict(
                title="Torque Feedback (lbf)",
                range=[
                    TORQUE_FEEDBACK_MIN_LBF,
                    TORQUE_FEEDBACK_MAX_LBF,
                ],
            ),
            aspectmode=aspect_mode,
        ),
        margin=dict(l=0, r=0, b=0, t=75),
    )

    plot_html = figure.to_html(
        full_html=False,
        include_plotlyjs=True,
        config={"responsive": True},
        div_id="cvtPlot",
    )

    all_data = {
        "rpm": points_rpm,
        "shift": points_shift,
        "torque": points_torque,
        "rpm_values": rpm_values,
        "shift_values": shift_values,
        "torque_values": torque_values,
    }

    constants = {
        "belt_width": BELT_WIDTH,
        "low_ratio": LOW_RATIO,
        "high_ratio": HIGH_RATIO,
        "primary_force_factor": PRIMARY_FORCE_FACTOR,
        "primary_ramp_reference_deg": PRIMARY_RAMP_REFERENCE_DEG,
        "torque_fit_a": TORQUE_FIT_A,
        "torque_fit_b": TORQUE_FIT_B,
        "torque_fit_c": TORQUE_FIT_C,
    }

    defaults = {
        "primary_spring_rate": PRIMARY_SPRING_RATE,
        "primary_pretension": PRIMARY_PRETENSION,
        "secondary_spring_rate": SECONDARY_SPRING_RATE,
        "secondary_pretension": SECONDARY_PRETENSION,
        "primary_weights": PRIMARY_WEIGHTS,
        "primary_weight_radius": PRIMARY_WEIGHT_RADIUS_IN,
        "primary_ramp_angle": PRIMARY_RAMP_ANGLE_DEG,
    }

    controls = "\n".join([
        model_input(
            "primarySpringRate",
            "Primary Spring Rate",
            PRIMARY_SPRING_RATE,
            0, 200, 1, "lbf/in",
        ),
        model_input(
            "primaryPretension",
            "Primary Pretension",
            PRIMARY_PRETENSION,
            0, 5, 0.01, "in",
        ),
        model_input(
            "secondarySpringRate",
            "Secondary Spring Rate",
            SECONDARY_SPRING_RATE,
            0, 100, 1, "lbf/in",
        ),
        model_input(
            "secondaryPretension",
            "Secondary Pretension",
            SECONDARY_PRETENSION,
            0, 5, 0.01, "in",
        ),
        model_input(
            "primaryWeights",
            "Primary Flyweight Total",
            PRIMARY_WEIGHTS,
            0, 5, 0.001, "lbf",
        ),
        model_input(
            "primaryWeightRadius",
            "Primary Weight Radius",
            PRIMARY_WEIGHT_RADIUS_IN,
            0, 5, 0.01, "in",
        ),
        model_input(
            "primaryRampAngle",
            "Primary Ramp Angle",
            PRIMARY_RAMP_ANGLE_DEG,
            1, 89, 0.1, "deg",
        ),
    ])

    selected = {
        "full": "",
        "slices": "",
        "cutaway": "",
    }
    selected[DEFAULT_VIEW_MODE] = "selected"

    sidebar_html = f"""
    <aside id="sidebar">
      <section class="sidebar-section">
        <h2>Model Variables</h2>
        <p class="small-note">
          Changes recalculate every visible point immediately.
          Normal values are the current Python values.
        </p>
        {controls}
        <button id="resetModelButton" class="full-button">
          Reset Model Values
        </button>
      </section>

      <section class="sidebar-section">
        <h2>Cube View</h2>

        <div class="model-row">
          <label for="viewMode">View Mode</label>
          <select id="viewMode">
            <option value="full" {selected["full"]}>Full Cube</option>
            <option value="slices" {selected["slices"]}>
              Cross-Section Slices
            </option>
            <option value="cutaway" {selected["cutaway"]}>
              Cutaway Cube
            </option>
          </select>
        </div>

        <p id="viewDescription" class="small-note"></p>

        <div class="slider-row">
          <label for="rpmCutSlider">RPM Cut / Slice</label>
          <input id="rpmCutSlider" type="range"
                 min="0" max="{len(rpm_values) - 1}"
                 value="{rpm_cut_index}" step="1">
          <span id="rpmCutValue"></span>
        </div>

        <div class="slider-row">
          <label for="shiftCutSlider">Shift Position Cut / Slice</label>
          <input id="shiftCutSlider" type="range"
                 min="0" max="{len(shift_values) - 1}"
                 value="{shift_cut_index}" step="1">
          <span id="shiftCutValue"></span>
        </div>

        <div class="slider-row">
          <label for="torqueCutSlider">Torque Feedback Cut / Slice</label>
          <input id="torqueCutSlider" type="range"
                 min="0" max="{len(torque_values) - 1}"
                 value="{torque_cut_index}" step="1">
          <span id="torqueCutValue"></span>
        </div>

        <p id="pointCount" class="small-note"></p>
      </section>

      <section class="sidebar-section">
        <h2>Data</h2>
        <button id="downloadCsvButton" class="full-button">
          Download Current Cube CSV
        </button>
        <p class="small-note">
          The downloaded CSV uses the model values currently shown
          in this sidebar.
        </p>
      </section>

      <section class="sidebar-section">
        <h2>Engine Torque Curve</h2>
        <p class="small-note">
          A constrained cubic least-squares fit is used:
        </p>
        <p class="equation-note">
          Torque = a × RPM³ + b × RPM² + c × RPM
        </p>
        <p class="small-note">
          It is forced through 0 RPM, 0 ft-lb. Values outside the
          measured 2400–3600 RPM dyno range are extrapolated.
        </p>
      </section>
    </aside>
    """

    javascript = r"""
    <script>
      const allData = __ALL_DATA_JSON__;
      const constants = __CONSTANTS_JSON__;
      const defaults = __DEFAULTS_JSON__;

      const plot = document.getElementById("cvtPlot");

      const viewMode = document.getElementById("viewMode");
      const rpmCutSlider = document.getElementById("rpmCutSlider");
      const shiftCutSlider = document.getElementById("shiftCutSlider");
      const torqueCutSlider = document.getElementById("torqueCutSlider");

      const rpmCutValue = document.getElementById("rpmCutValue");
      const shiftCutValue = document.getElementById("shiftCutValue");
      const torqueCutValue = document.getElementById("torqueCutValue");
      const viewDescription = document.getElementById("viewDescription");
      const pointCount = document.getElementById("pointCount");

      const inputIds = [
        "primarySpringRate",
        "primaryPretension",
        "secondarySpringRate",
        "secondaryPretension",
        "primaryWeights",
        "primaryWeightRadius",
        "primaryRampAngle"
      ];

      let currentNetForce = [];
      let currentCustomData = [];

      function numberFromInput(id, fallback, minimum) {
        const value = Number(document.getElementById(id).value);

        if (!Number.isFinite(value)) {
          return fallback;
        }

        return Math.max(minimum, value);
      }

      function currentModel() {
        return {
          primarySpringRate: numberFromInput(
            "primarySpringRate",
            defaults.primary_spring_rate,
            0
          ),
          primaryPretension: numberFromInput(
            "primaryPretension",
            defaults.primary_pretension,
            0
          ),
          secondarySpringRate: numberFromInput(
            "secondarySpringRate",
            defaults.secondary_spring_rate,
            0
          ),
          secondaryPretension: numberFromInput(
            "secondaryPretension",
            defaults.secondary_pretension,
            0
          ),
          primaryWeights: numberFromInput(
            "primaryWeights",
            defaults.primary_weights,
            0
          ),
          primaryWeightRadius: numberFromInput(
            "primaryWeightRadius",
            defaults.primary_weight_radius,
            0
          ),
          primaryRampAngle: numberFromInput(
            "primaryRampAngle",
            defaults.primary_ramp_angle,
            0.1
          )
        };
      }

      function engineTorqueAtRPM(rpm) {
        const torque = (
          constants.torque_fit_a * rpm * rpm * rpm
          + constants.torque_fit_b * rpm * rpm
          + constants.torque_fit_c * rpm
        );

        return Math.max(0, torque);
      }

      function ratioAtShift(shift) {
        let fraction = shift / constants.belt_width;
        fraction = Math.max(0, Math.min(1, fraction));

        return constants.low_ratio + fraction * (
          constants.high_ratio - constants.low_ratio
        );
      }

      function rampMultiplier(angleDeg) {
        const refRadians = (
          constants.primary_ramp_reference_deg * Math.PI / 180
        );
        const angleRadians = angleDeg * Math.PI / 180;
        const refSine = Math.sin(refRadians);

        if (Math.abs(refSine) < 1e-12) {
          return 1;
        }

        return Math.sin(angleRadians) / refSine;
      }

      function calculatePoint(rpm, shift, torqueFeedback, model) {
        const ratio = ratioAtShift(shift);
        const engineTorque = engineTorqueAtRPM(rpm);

        const flyweightForce = (
          rpm
          * model.primaryWeights
          * model.primaryWeightRadius
          * constants.primary_force_factor
          * rampMultiplier(model.primaryRampAngle)
        );

        const primarySpringForce = (
          model.primarySpringRate
          * (model.primaryPretension - shift)
        );

        const primaryForce = flyweightForce - primarySpringForce;

        const secondarySpring = (
          model.secondarySpringRate
          * (model.secondaryPretension + shift)
        );

        const secondaryTotal = secondarySpring + torqueFeedback;
        const netForce = primaryForce - secondaryTotal;

        return {
          ratio,
          engineTorque,
          primaryForce,
          secondarySpring,
          secondaryTotal,
          netForce
        };
      }

      function recalculateAllPoints() {
        const model = currentModel();

        currentNetForce = new Array(allData.rpm.length);
        currentCustomData = new Array(allData.rpm.length);

        for (let i = 0; i < allData.rpm.length; i++) {
          const result = calculatePoint(
            allData.rpm[i],
            allData.shift[i],
            allData.torque[i],
            model
          );

          currentNetForce[i] = result.netForce;

          currentCustomData[i] = [
            result.ratio,
            result.engineTorque,
            result.primaryForce,
            result.secondarySpring,
            result.secondaryTotal,
            result.netForce
          ];
        }
      }

      function selectedCuts() {
        return {
          rpm: allData.rpm_values[Number(rpmCutSlider.value)],
          shift: allData.shift_values[Number(shiftCutSlider.value)],
          torque: allData.torque_values[Number(torqueCutSlider.value)]
        };
      }

      function updateLabels() {
        const cuts = selectedCuts();

        rpmCutValue.textContent = cuts.rpm.toFixed(0) + " RPM";
        shiftCutValue.textContent = cuts.shift.toFixed(3) + " in";
        torqueCutValue.textContent = cuts.torque.toFixed(1) + " lbf";

        if (viewMode.value === "full") {
          viewDescription.textContent =
            "Showing every RPM, shift-position, and torque-feedback combination.";
        } else if (viewMode.value === "slices") {
          viewDescription.textContent =
            "Showing three internal cut planes at the selected RPM, shift position, and torque feedback.";
        } else {
          viewDescription.textContent =
            "Showing only points at or below the selected RPM, shift position, and torque feedback.";
        }
      }

      function nearlyEqual(a, b) {
        return Math.abs(a - b) < 0.000001;
      }

      function updateVisiblePoints() {
        const cuts = selectedCuts();
        const mode = viewMode.value;

        const x = [];
        const y = [];
        const z = [];
        const colors = [];
        const customData = [];

        for (let i = 0; i < allData.rpm.length; i++) {
          const rpm = allData.rpm[i];
          const shift = allData.shift[i];
          const torque = allData.torque[i];

          let visible = true;

          if (mode === "slices") {
            visible = (
              nearlyEqual(rpm, cuts.rpm)
              || nearlyEqual(shift, cuts.shift)
              || nearlyEqual(torque, cuts.torque)
            );
          } else if (mode === "cutaway") {
            visible = (
              rpm <= cuts.rpm
              && shift <= cuts.shift
              && torque <= cuts.torque
            );
          }

          if (visible) {
            x.push(rpm);
            y.push(shift);
            z.push(torque);
            colors.push(currentNetForce[i]);
            customData.push(currentCustomData[i]);
          }
        }

        Plotly.restyle(
          plot,
          {
            x: [x],
            y: [y],
            z: [z],
            "marker.color": [colors],
            customdata: [customData]
          },
          [0]
        );

        pointCount.textContent =
          "Visible points: "
          + x.length.toLocaleString()
          + " of "
          + allData.rpm.length.toLocaleString();
      }

      function updateEverything() {
        recalculateAllPoints();
        updateLabels();
        updateVisiblePoints();
      }

      function resetModel() {
        document.getElementById("primarySpringRate").value =
          defaults.primary_spring_rate;
        document.getElementById("primaryPretension").value =
          defaults.primary_pretension;
        document.getElementById("secondarySpringRate").value =
          defaults.secondary_spring_rate;
        document.getElementById("secondaryPretension").value =
          defaults.secondary_pretension;
        document.getElementById("primaryWeights").value =
          defaults.primary_weights;
        document.getElementById("primaryWeightRadius").value =
          defaults.primary_weight_radius;
        document.getElementById("primaryRampAngle").value =
          defaults.primary_ramp_angle;

        updateEverything();
      }

      function csvValue(value) {
        const text = String(value);

        if (
          text.includes(",")
          || text.includes('"')
          || text.includes("\n")
        ) {
          return '"' + text.replaceAll('"', '""') + '"';
        }

        return text;
      }

      function downloadCurrentCsv() {
        const header = [
          "rpm",
          "shift_position_in",
          "torque_feedback_lbf",
          "cvt_ratio",
          "engine_torque_ftlb",
          "primary_force_lbf",
          "secondary_spring_lbf",
          "secondary_total_lbf",
          "net_shift_force_lbf"
        ];

        const lines = [header.join(",")];

        for (let i = 0; i < allData.rpm.length; i++) {
          const calculated = currentCustomData[i];

          const row = [
            allData.rpm[i],
            allData.shift[i],
            allData.torque[i],
            calculated[0],
            calculated[1],
            calculated[2],
            calculated[3],
            calculated[4],
            calculated[5]
          ];

          lines.push(row.map(csvValue).join(","));
        }

        const blob = new Blob(
          [lines.join("\n")],
          {type: "text/csv;charset=utf-8;"}
        );

        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");

        link.href = url;
        link.download = "cvt_torque_feedback_current_model.csv";

        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);

        URL.revokeObjectURL(url);
      }

      inputIds.forEach(function(id) {
        const input = document.getElementById(id);
        input.addEventListener("input", updateEverything);
        input.addEventListener("change", updateEverything);
      });

      viewMode.addEventListener("change", updateEverything);
      rpmCutSlider.addEventListener("input", updateEverything);
      shiftCutSlider.addEventListener("input", updateEverything);
      torqueCutSlider.addEventListener("input", updateEverything);

      document.getElementById("resetModelButton").addEventListener(
        "click",
        resetModel
      );

      document.getElementById("downloadCsvButton").addEventListener(
        "click",
        downloadCurrentCsv
      );

      setTimeout(updateEverything, 100);
    </script>
    """

    javascript = javascript.replace(
        "__ALL_DATA_JSON__",
        json.dumps(all_data),
    )
    javascript = javascript.replace(
        "__CONSTANTS_JSON__",
        json.dumps(constants),
    )
    javascript = javascript.replace(
        "__DEFAULTS_JSON__",
        json.dumps(defaults),
    )

    page = f"""
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="utf-8">
      <title>CVT Torque Feedback Parameter Sweep</title>
      <style>
        * {{
          box-sizing: border-box;
        }}

        body {{
          margin: 0;
          background: #ffffff;
          color: #1f2937;
          font-family: Arial, Helvetica, sans-serif;
        }}

        #app {{
          min-height: 100vh;
          display: grid;
          grid-template-columns: 340px minmax(0, 1fr);
        }}

        #sidebar {{
          padding: 18px;
          overflow-y: auto;
          background: #f8fafc;
          border-right: 1px solid #cbd5e1;
        }}

        #plotArea {{
          min-width: 0;
          padding: 8px;
        }}

        .sidebar-section {{
          margin-bottom: 18px;
          padding: 14px;
          background: #ffffff;
          border: 1px solid #dbe3ec;
          border-radius: 10px;
        }}

        .sidebar-section h2 {{
          margin: 0 0 12px 0;
          font-size: 18px;
        }}

        .model-row,
        .slider-row {{
          margin-bottom: 13px;
        }}

        .model-row label,
        .slider-row label {{
          display: block;
          margin-bottom: 5px;
          font-size: 13px;
          font-weight: 600;
        }}

        .input-group {{
          display: grid;
          grid-template-columns: 1fr 70px;
          gap: 7px;
          align-items: center;
        }}

        .input-group span {{
          color: #475569;
          font-size: 12px;
        }}

        input[type="number"],
        select {{
          width: 100%;
          padding: 7px;
          border: 1px solid #94a3b8;
          border-radius: 6px;
          background: #ffffff;
        }}

        input[type="range"] {{
          width: 100%;
        }}

        .slider-row span {{
          display: block;
          margin-top: 3px;
          color: #475569;
          font-size: 12px;
          font-variant-numeric: tabular-nums;
        }}

        .full-button {{
          width: 100%;
          padding: 9px;
          border: 1px solid #64748b;
          border-radius: 6px;
          background: #ffffff;
          color: #0f172a;
          font-weight: 600;
          cursor: pointer;
        }}

        .full-button:hover {{
          background: #e2e8f0;
        }}

        .small-note {{
          margin: 8px 0 0 0;
          color: #475569;
          font-size: 12px;
          line-height: 1.45;
        }}

        .equation-note {{
          margin: 8px 0;
          padding: 8px;
          border-radius: 5px;
          background: #f1f5f9;
          font-family: Consolas, monospace;
          font-size: 12px;
        }}

        @media (max-width: 900px) {{
          #app {{
            grid-template-columns: 1fr;
          }}

          #sidebar {{
            border-right: none;
            border-bottom: 1px solid #cbd5e1;
          }}
        }}
      </style>
    </head>
    <body>
      <div id="app">
        {sidebar_html}
        <main id="plotArea">
          {plot_html}
        </main>
      </div>
      {javascript}
    </body>
    </html>
    """

    output_path = Path(OUTPUT_FILE_NAME).resolve()
    output_path.write_text(page, encoding="utf-8")

    if SAVE_DEFAULT_CSV:
        csv_path = write_csv(default_rows)
        print(f"Default calculation CSV saved: {csv_path}")

    print("\nGraph created successfully.")
    print(
        "Torque fit: "
        f"{TORQUE_FIT_A:.12e}*RPM^3 + "
        f"{TORQUE_FIT_B:.12e}*RPM^2 + "
        f"{TORQUE_FIT_C:.12e}*RPM"
    )
    print(f"RPM values: {len(rpm_values)}")
    print(f"Shift positions: {len(shift_values)}")
    print(f"Torque-feedback layers: {len(torque_values)}")
    print(f"Total cube points: {len(default_rows)}")
    print(f"Interactive graph saved: {output_path}")

    webbrowser.open(output_path.as_uri())


if __name__ == "__main__":
    make_graph()
