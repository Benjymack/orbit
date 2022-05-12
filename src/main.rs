use std::f64::consts::PI;
use nalgebra::Vector2;
use macroquad::prelude::*;
use std::collections::BTreeSet;

// Types
type Displacement = Vector2<f64>;
type Velocity = Vector2<f64>;
type Force = Vector2<f64>;

// Constants
const DT: f64 = 20.; // s
const G: f64 = 6.67384e-11; // m^3/kg/s^2
const PLANET_DENSITY: f64 = 5520.;

fn sphere_volume(radius: f64) -> f64 {
    4./3. * PI * radius.powi(3)
}
fn sphere_radius(volume: f64) -> f64 {
    (volume * 3./4. / PI).cbrt()
}

fn gravity_force(r1: &Displacement, r2: &Displacement, m1: f64, m2: f64) -> Force {
    // NOTE: Calculates the force on r1! (negate to get force on r2)
    let displacement = r1 - r2;
    let force_magnitude = G * m1 * m2 / displacement.magnitude_squared();
    force_magnitude * -displacement.normalize()
}

fn better_draw_circle(x: f64, y: f64, r: f64, color: Color) {
    draw_poly(x as f32, y as f32, 100, r as f32, 0., color);
}

// Body
#[derive(Copy, Clone)]
struct Body {
    displacement: Displacement,
    velocity: Velocity,
    mass: f64,
    radius: f64,
    color: Color,
}

impl Body {
    fn update(&mut self, force: &Force) {
        let acceleration = force / self.mass;
        self.velocity += acceleration * DT;
        self.displacement += self.velocity * DT;
    }

    fn new_planet_radius(displacement: Displacement, velocity: Velocity, radius: f64) -> Body {
        Body {
            displacement,
            velocity,
            mass: sphere_volume(radius) * PLANET_DENSITY,
            radius,
            color: PURPLE,
        }
    }

    fn new_planet_mass(displacement: Displacement, velocity: Velocity, mass: f64) -> Body {
        Body {
            displacement,
            velocity,
            mass,
            radius: sphere_radius(mass / PLANET_DENSITY),
            color: PURPLE,
        }
    }
}

enum BodyCreation {
    None,
    Position(Body), // x and y
    Velocity(Body), // x y and z
}

// Functions
fn draw_body(body: &Body, origin_x: f64, origin_y: f64, scale_factor: f64) {
    let canvas_x = origin_x + scale_factor * body.displacement.x;
    let canvas_y = origin_y + scale_factor * body.displacement.y;
    let canvas_radius = scale_factor * body.radius;

    better_draw_circle(canvas_x, canvas_y, canvas_radius, body.color);
}

fn get_real_mouse_pos(origin_x: f64, origin_y: f64, scale_factor: f64) -> Displacement {
    let real_x = (mouse_position().0 as f64 - origin_x) / scale_factor;
    let real_y = (mouse_position().1 as f64 - origin_y) / scale_factor;

    Displacement::new(real_x, real_y)
}

fn draw_arrow(x1: f32, y1: f32, x2: f32, y2: f32, color: Color) {
    const THICKNESS: f32 = 1.;
    const ARROWHEAD_ANGLE: f32 = (PI / 8.) as f32;
    const ARROWHEAD_LENGTH: f32 = 10.;

    draw_line(x1, y1, x2, y2, THICKNESS, color);

    let initial_angle = (y2 - y1).atan2(x2 - x1);
    let angle1 = initial_angle + ARROWHEAD_ANGLE;
    let angle2 = initial_angle - ARROWHEAD_ANGLE;

    let arrowhead1_x = x2 - ARROWHEAD_LENGTH * angle1.cos();
    let arrowhead1_y = y2 - ARROWHEAD_LENGTH * angle1.sin();
    let arrowhead2_x = x2 - ARROWHEAD_LENGTH * angle2.cos();
    let arrowhead2_y = y2 - ARROWHEAD_LENGTH * angle2.sin();

    draw_line(x2, y2, arrowhead1_x, arrowhead1_y, THICKNESS, color);
    draw_line(x2, y2, arrowhead2_x, arrowhead2_y, THICKNESS, color);
}

fn draw_velocity(body: &Body, origin_x: f64, origin_y: f64, scale_factor: f64) {
    const ARBITRARY_CONSTANT: f64 = 500.;

    if body.velocity == Velocity::zeros() {
        return;
    }

    let x1 = origin_x + scale_factor * body.displacement.x;
    let y1 = origin_y + scale_factor * body.displacement.y;
    let x2 = x1 + scale_factor * body.velocity.x * ARBITRARY_CONSTANT;
    let y2 = y1 + scale_factor * body.velocity.y * ARBITRARY_CONSTANT;

    draw_arrow(x1 as f32, y1 as f32, x2 as f32, y2 as f32, GREEN);
}

fn calculate_trajectory(bodies: &Vec<Body>, target_body: &Body) -> (Vec<Displacement>, bool) {
    let mut trajectory: Vec<Displacement> = Vec::new();

    let mut temp_bodies = vec![*target_body];
    temp_bodies.extend(bodies);

    let mut does_collide = false;

    for _ in 0..1000 {
        let removed = update(&mut temp_bodies);

        if removed.contains(&0) {
            does_collide = true;
            break;
        }

        trajectory.push(temp_bodies[0].displacement);
    }

    (trajectory, does_collide)
}

fn draw(bodies: &Vec<Body>, is_paused: &bool, current_body: &BodyCreation, origin_x: f64, origin_y: f64, scale_factor: f64, explosion_texture: &Texture2D) {
    // Draw planets
    for body in bodies {
        draw_body(body, origin_x, origin_y, scale_factor);
        draw_velocity(body, origin_x, origin_y, scale_factor);
    }
    match current_body {
        BodyCreation::Position(body) => {
            draw_body(body, origin_x, origin_y, scale_factor);
        },
        BodyCreation::Velocity(body) => {
            draw_body(body, origin_x, origin_y, scale_factor);
            draw_arrow(
                (origin_x + body.displacement.x * scale_factor) as f32,
                (origin_y + body.displacement.y * scale_factor) as f32,
                mouse_position().0,
                mouse_position().1, WHITE);

            // Draw the trajectory
            const TRAJECTORY_THICKNESS: f32 = 1.;
            const TRAJECTORY_COLOR: Color = YELLOW;
            let (trajectory, does_collide) = calculate_trajectory(bodies, body);
            if trajectory.len() > 0 {
                for i in 0..(trajectory.len() - 1) {
                    let x1 = (origin_x + trajectory[i].x * scale_factor) as f32;
                    let y1 = (origin_y + trajectory[i].y * scale_factor) as f32;
                    let x2 = (origin_x + trajectory[i + 1].x * scale_factor) as f32;
                    let y2 = (origin_y + trajectory[i + 1].y * scale_factor) as f32;
                    draw_line(x1, y1, x2, y2, TRAJECTORY_THICKNESS, TRAJECTORY_COLOR);
                }
                if does_collide {
                    let x = (origin_x + trajectory[trajectory.len() - 1].x * scale_factor) as f32;
                    let y = (origin_y + trajectory[trajectory.len() - 1].y * scale_factor) as f32;
                    draw_texture(*explosion_texture, x - explosion_texture.width() / 2., y - explosion_texture.height() / 2., WHITE);
                }
            }
        },
        BodyCreation::None => {},
    };

    // Draw running/paused state
    let (text, color) = if *is_paused {
        ("Paused", RED)
    } else {
        ("Running", GREEN)
    };
    const TEXT_X: f32 = 10.;
    const TEXT_Y: f32 = 10.;
    const TEXT_SIZE: f32 = 14.;
    draw_text(text, TEXT_X, TEXT_Y+TEXT_SIZE, TEXT_SIZE, color);

    // Draw time speed
    draw_text(&*format!("{:.3}x speed", DT * (get_fps() as f64)), TEXT_X, TEXT_Y+2.*TEXT_SIZE, TEXT_SIZE, WHITE);

    // Draw scale
    const LINE_X: f32 = 10.;
    const LINE_LENGTH: f32 = 100.;
    const LINE_Y: f32 = 10.;
    const SIDE_LENGTH: f32 = 5.;
    const LINE_THICKNESS: f32 = 1.;

    draw_line(LINE_X, screen_height()-LINE_Y-SIDE_LENGTH, LINE_X, screen_height()-LINE_Y+SIDE_LENGTH, LINE_THICKNESS, WHITE);
    draw_line(LINE_X+ LINE_LENGTH, screen_height()-LINE_Y-SIDE_LENGTH, LINE_X+ LINE_LENGTH, screen_height()-LINE_Y+SIDE_LENGTH, LINE_THICKNESS, WHITE);
    draw_line(LINE_X, screen_height()-LINE_Y, LINE_X+ LINE_LENGTH, screen_height()-LINE_Y, LINE_THICKNESS, WHITE);

    // Draw the scale text
    let scale_text = format!("{:.3} m", SIDE_LENGTH as f64 / scale_factor);
    draw_text(&*scale_text, LINE_X+5., screen_height()-LINE_Y-10., 12., WHITE);
}

fn update(bodies: &mut Vec<Body>) -> BTreeSet<usize> {
    // Calculate the forces
    let mut forces = vec![Force::zeros(); bodies.len()];

    for i in 0..bodies.len() {
        let mut net_force = Force::zeros();
        for j in 0..bodies.len() {
            if i == j {
                continue;
            }
            net_force += gravity_force(&bodies[i].displacement, &bodies[j].displacement,
                                       bodies[i].mass, bodies[j].mass);
        }
        forces[i] = net_force;
    }

    // Apply the forces
    for i in 0..bodies.len() {
        bodies[i].update(&forces[i]);
    }

    // Collisions
    let mut to_remove = BTreeSet::new();

    for i in 1..bodies.len() {
        for j in 0..i {
            if (bodies[i].displacement - bodies[j].displacement).magnitude() <= bodies[i].radius + bodies[j].radius {
                // They collide
                let new_mass = bodies[i].mass + bodies[j].mass;
                let new_velocity = (bodies[i].velocity * bodies[i].mass + bodies[j].velocity * bodies[j].mass) / new_mass;
                let new_displacement = if bodies[i].mass > bodies[j].mass {
                    bodies[i].displacement
                } else {
                    bodies[j].displacement
                };

                to_remove.insert(i);
                to_remove.insert(j);

                bodies.push(Body::new_planet_mass(new_displacement, new_velocity, new_mass));
            }
        }
    }

    // Remove the bodies that have collided
    for idx in to_remove.iter().rev() {
        bodies.remove(*idx);
    }

    to_remove
}

fn process_input(bodies: &mut Vec<Body>, is_paused: &mut bool, current_body: &mut BodyCreation, origin_x: f64, origin_y: f64, scale_factor: f64, autoscale: &mut bool, manual_scale_factor: &mut f64, manual_origin_x: &mut f64, manual_origin_y: &mut f64) {
    // Pause
    if is_key_pressed(KeyCode::Space) {
        *is_paused = !*is_paused;
    }

    // Remove all bodies
    if is_key_pressed(KeyCode::C) {
        bodies.clear();
    }

    // Remove a particular body
    if is_mouse_button_pressed(MouseButton::Right) {
        for (i, body) in bodies.iter().enumerate() {
            if (get_real_mouse_pos(origin_x, origin_y, scale_factor) - body.displacement).magnitude() < body.radius {
                bodies.remove(i);
                break;
            }
        }
    }

    // Create a new body
    if is_key_pressed(KeyCode::Escape) {
        *current_body = BodyCreation::None;
    }

    if is_key_down(KeyCode::LeftControl) {
        if is_key_down(KeyCode::Equal) {
            if *autoscale {
                *autoscale = false;
                *manual_scale_factor = scale_factor;
                *manual_origin_x = origin_x;
                *manual_origin_y = origin_y;
            }
            *manual_scale_factor *= 1.1;
        } else if is_key_down(KeyCode::Minus) {
            if *autoscale {
                *autoscale = false;
                *manual_scale_factor = scale_factor;
                *manual_origin_x = origin_x;
                *manual_origin_y = origin_y;
            }
            *manual_scale_factor /= 1.1;
        }
    }

    if is_key_down(KeyCode::Left) {
        if *autoscale {
            *autoscale = false;
            *manual_scale_factor = scale_factor;
            *manual_origin_x = origin_x;
            *manual_origin_y = origin_y;
        }
        *manual_origin_x += 10.;
    } else if is_key_down(KeyCode::Right) {
        if *autoscale {
            *autoscale = false;
            *manual_scale_factor = scale_factor;
            *manual_origin_x = origin_x;
            *manual_origin_y = origin_y;
        }
        *manual_origin_x -= 10.;
    }
    if is_key_down(KeyCode::Up) {
        if *autoscale {
            *autoscale = false;
            *manual_scale_factor = scale_factor;
            *manual_origin_x = origin_x;
            *manual_origin_y = origin_y;
        }
        *manual_origin_y += 10.;
    } else if is_key_down(KeyCode::Down) {
        if *autoscale {
            *autoscale = false;
            *manual_scale_factor = scale_factor;
            *manual_origin_x = origin_x;
            *manual_origin_y = origin_y;
        }
        *manual_origin_y -= 10.;
    }

    match current_body {
        BodyCreation::None => {
            if is_mouse_button_pressed(MouseButton::Left) {
                *current_body = BodyCreation::Position(Body::new_planet_radius(
                    get_real_mouse_pos(origin_x, origin_y, scale_factor),
                    Velocity::zeros(),
                    0.));
            }
        }
        BodyCreation::Position(body) => {
            if is_mouse_button_down(MouseButton::Left) {
                let radius = (get_real_mouse_pos(origin_x, origin_y, scale_factor) - body.displacement).magnitude();
                body.radius = radius;
            } else {
                if body.radius == 0. {
                    *current_body = BodyCreation::None;
                } else {
                    *current_body = BodyCreation::Velocity(Body::new_planet_radius(body.displacement, body.velocity, body.radius));
                }
            }
        }
        BodyCreation::Velocity(body) => {
            // Set the velocity based on the mouse position
            let velocity_displacement = get_real_mouse_pos(origin_x, origin_y, scale_factor) - body.displacement;
            let velocity = velocity_displacement / 1000.; // Have to convert the displacement to a velocity with an arbitrary constant, representing the time to travel that distance

            body.velocity = velocity;

            if is_key_pressed(KeyCode::F) {
                // Set the velocity to zero
                body.velocity = Velocity::zeros();
                bodies.push(*body);
                *current_body = BodyCreation::None;
            } else if is_mouse_button_pressed(MouseButton::Left) {
                // Finalise the body
                bodies.push(*body);
                *current_body = BodyCreation::None;
            }
        }
    };
}

#[macroquad::main("Orbit")]
async fn main() {
    let mut bodies: Vec<Body> = Vec::new();

    let mut is_paused = false;
    let mut current_body = BodyCreation::None;

    let mut autoscale = true;

    let mut manual_scale_factor = 1e-4;
    let mut manual_origin_x = 0.;
    let mut manual_origin_y = 0.;

    let explosion_texture = Texture2D::from_image(&Image::from_file_with_format(include_bytes!("../explosion.png"), Some(ImageFormat::Png)));

    loop {
        // Determine scale
        let (scale_factor, origin_x, origin_y) = if autoscale {
            // Autoscale will automatically grow to display all of the planets
            let mut min_x = 0_f64;
            let mut min_y = 0_f64;
            let mut max_x = 0_f64;
            let mut max_y = 0_f64;

            if bodies.is_empty() {
                min_x = -1000.;
                min_y = min_x;
                max_x = -min_x;
                max_y = -min_x;
            } else {
                for body in &bodies {
                    min_x = min_x.min(body.displacement.x - body.radius);
                    max_x = max_x.max(body.displacement.x + body.radius);

                    min_y = min_y.min(body.displacement.y - body.radius);
                    max_y = max_y.max(body.displacement.y + body.radius);
                }
            }

            const EXTRA: f64 = 1.1;

            let half_width = max_x.max(-min_x) * EXTRA;
            let half_height = max_y.max(-min_y) * EXTRA;

            // Calculate the origin (top-left corner) and the scale factor
            let required_aspect_ratio = half_width / half_height;
            let actual_aspect_ratio = screen_width() / screen_height();

            // The scale factor is pixels per meter
            // The origin x and y are the position (in pixels) of the origin point
            let origin_x = screen_width() as f64 / 2.;
            let origin_y = screen_height() as f64 / 2.;


            let scale_factor = if required_aspect_ratio > actual_aspect_ratio as f64 {
                // We will have extra space at the top and bottom
                screen_width() as f64 / half_width
            } else {
                // We will have extra space on the left and right
                screen_height() as f64 / half_height
            } / 2.;

            (scale_factor, origin_x, origin_y)
        } else {
            // Manual scale will not change automatically

            // let origin_x = screen_width() as f64 / 2.;
            // let origin_y = screen_height() as f64 / 2.;

            (manual_scale_factor, manual_origin_x, manual_origin_y)
        };


        process_input(&mut bodies, &mut is_paused, &mut current_body, origin_x, origin_y, scale_factor, &mut autoscale, &mut manual_scale_factor, &mut manual_origin_x, &mut manual_origin_y);
        if !is_paused {
            update(&mut bodies);
        }
        draw(&bodies, &is_paused, &current_body, origin_x, origin_y, scale_factor, &explosion_texture);
        next_frame().await
    }
}
