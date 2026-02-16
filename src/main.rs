use macroquad::prelude::*;
use std::f32::consts::TAU;

#[macroquad::main("Braideroids : Asteroids = Braid;")]
async fn main() {
    let mut asteroid1 = Asteroid::new(Body::new(250.0, 250.0, -20.0, 10.0, -2.0, 0.5), 3, 120.0);
    let mut asteroid2 = Asteroid::new(Body::new(250.0, 250.0, 30.0, -30.0, 4.0, -1.0), 3, 50.0);
    let mut asteroid3 = Asteroid::new(Body::new(250.0, 250.0, 40.0, -50.0, 5.0, -0.3), 3, 50.0);
    let mut asteroid4 = Asteroid::new(Body::new(250.0, 250.0, -50.0, -10.0, 3.0, 0.5), 5, 100.0);

    let mut game_last_tick = get_time() as f32;
    loop {
        let current_tick = get_time() as f32;
        let dt = current_tick - game_last_tick;
        game_last_tick = current_tick;

        asteroid1.body.update(dt);
        asteroid2.body.update(dt);
        asteroid3.body.update(dt);
        asteroid4.body.update(dt);

        reflect(&mut asteroid1.body);
        reflect(&mut asteroid2.body);
        reflect(&mut asteroid3.body);
        reflect(&mut asteroid4.body);

        asteroid1.draw();
        asteroid2.draw();
        asteroid3.draw();
        asteroid4.draw();

        highlight_collision(&asteroid1, &asteroid2);
        highlight_collision(&asteroid1, &asteroid3);
        highlight_collision(&asteroid1, &asteroid4);

        highlight_collision(&asteroid2, &asteroid3);
        highlight_collision(&asteroid2, &asteroid4);

        highlight_collision(&asteroid3, &asteroid4);

        next_frame().await;
    }
}

fn draw_vertices(object: &impl Shape) {
    let vertices1 = object.shape();
    vertices1.iter().for_each(|vertice| {
        draw_circle(vertice.x, vertice.y, 5.0, BLUE);
    });
}

fn highlight_collision(asteroid1: &Asteroid, asteroid2: &Asteroid) {
    let collided = collision(asteroid1, asteroid2);
    if let Some(collision_info) = collided {
        /* draw_poly(
            asteroid1.body.lin_pos.x,
            asteroid1.body.lin_pos.y,
            asteroid1.sides,
            asteroid1.size,
            asteroid1.body.ang_pos.to_degrees(),
            GRAY,
        );
        draw_poly(
            asteroid2.body.lin_pos.x,
            asteroid2.body.lin_pos.y,
            asteroid2.sides,
            asteroid2.size,
            asteroid2.body.ang_pos.to_degrees(),
            GRAY,
        ); */
        collision_info.which_object_which_vertice[0]
            .iter()
            .for_each(|vertice_id| {
                let vertice = asteroid1.shape()[*vertice_id];
                draw_circle(vertice.x, vertice.y, 4.0, RED);

                let vertice_p = asteroid1.shape()
                    [(*vertice_id + asteroid1.shape().len() - 1) % asteroid1.shape().len()];
                let vertice_n = asteroid1.shape()[(*vertice_id + 1) % asteroid1.shape().len()];
                draw_line(vertice.x, vertice.y, vertice_p.x, vertice_p.y, 3.0, RED);
                draw_line(vertice.x, vertice.y, vertice_n.x, vertice_n.y, 3.0, RED);
            });
        collision_info.which_object_which_vertice[1]
            .iter()
            .for_each(|vertice_id| {
                let vertice = asteroid2.shape()[*vertice_id];
                draw_circle(vertice.x, vertice.y, 4.0, RED);

                let vertice_p = asteroid2.shape()
                    [(*vertice_id + asteroid2.shape().len() - 1) % asteroid2.shape().len()];
                let vertice_n = asteroid2.shape()[(*vertice_id + 1) % asteroid2.shape().len()];
                draw_line(vertice.x, vertice.y, vertice_p.x, vertice_p.y, 3.0, RED);
                draw_line(vertice.x, vertice.y, vertice_n.x, vertice_n.y, 3.0, RED);
            });
    }
}

trait Shape {
    fn shape(&self) -> Vec<Vec2>;
}

impl Shape for Asteroid {
    fn shape(&self) -> Vec<Vec2> {
        let center = self.body.lin_pos;
        let radius = self.size;
        let rotation = self.body.ang_pos;

        let mut vertices = vec![];
        let theta = TAU / self.sides as f32;
        for it in 0..self.sides {
            vertices.push(
                center
                    + radius
                        * Vec2::new(
                            (it as f32 * theta + rotation).cos(),
                            (it as f32 * theta + rotation).sin(),
                        ),
            );
        }
        return vertices;
    }
}

struct CollisionInfo {
    which_object_which_vertice: [Vec<usize>; 2],
}

fn collision(object1: &impl Shape, object2: &impl Shape) -> Option<CollisionInfo> {
    let mut collided = false;

    let mut collision_info = CollisionInfo {
        which_object_which_vertice: [Vec::new(), Vec::new()],
    };

    for (it, vertice) in object1.shape().iter().enumerate() {
        if inside(&vertice, object2) {
            collision_info.which_object_which_vertice[0].push(it);
            collided = true;
        }
    }

    for (it, vertice) in object2.shape().iter().enumerate() {
        let is_inside = inside(&vertice, object1);
        if is_inside {
            collision_info.which_object_which_vertice[1].push(it);
            collided = true;
        }
    }

    if collided {
        return Some(collision_info);
    }
    None
}

fn rotate(vector: &Vec2, theta: f32) -> Vec2 {
    Mat2 {
        x_axis: Vec2 {
            x: theta.cos(),
            y: -1.0 * theta.sin(),
        },
        y_axis: Vec2 {
            x: theta.sin(),
            y: theta.cos(),
        },
    } * *vector
}

fn inside(vertice: &Vec2, object: &impl Shape) -> bool {
    let vertices = object.shape();
    let mut surface_perps = vec![];

    for it in 0..vertices.len() {
        let b = vertices[(it + 1) % vertices.len()];
        let a = vertices[it % vertices.len()];
        let normal = b - a;
        let orthogonal = rotate(&normal, TAU / 4.0);

        surface_perps.push(orthogonal);
    }

    for axis in surface_perps {
        let projections: Vec<f32> = vertices.iter().map(|vertice| vertice.dot(axis)).collect();
        let min_proj = (&projections)
            .into_iter()
            .min_by(|a, b| a.total_cmp(b))
            .unwrap();
        let max_proj = (&projections)
            .into_iter()
            .max_by(|a, b| a.total_cmp(b))
            .unwrap();
        let vertice_proj = vertice.dot(axis);

        if vertice_proj < *min_proj {
            return false;
        }

        if vertice_proj > *max_proj {
            return false;
        }
    }

    return true;
}

fn reflect(body: &mut Body) {
    if body.lin_pos.x >= screen_width() * 0.75 {
        body.lin_pos.x = screen_width() * 0.75;
        body.lin_vel = Vec2::new(-1.0 * body.lin_vel.x, body.lin_vel.y)
    }
    if body.lin_pos.x <= screen_width() * 0.25 {
        body.lin_pos.x = screen_width() * 0.25;
        body.lin_vel = Vec2::new(-1.0 * body.lin_vel.x, body.lin_vel.y)
    }
    if body.lin_pos.y >= screen_height() * 0.75 {
        body.lin_pos.y = screen_height() * 0.75;
        body.lin_vel = Vec2::new(body.lin_vel.x, -1.0 * body.lin_vel.y)
    }
    if body.lin_pos.y <= screen_height() * 0.25 {
        body.lin_pos.y = screen_height() * 0.25;
        body.lin_vel = Vec2::new(body.lin_vel.x, -1.0 * body.lin_vel.y)
    }
}

#[derive(Clone)]
struct Body {
    lin_pos: Vec2,
    lin_vel: Vec2,
    lin_acc: Vec2,

    ang_pos: f32,
    ang_vel: f32,
    ang_acc: f32,
}

impl Body {
    fn new(x0: f32, y0: f32, vx0: f32, vy0: f32, t0: f32, o0: f32) -> Self {
        Self {
            lin_pos: Vec2 { x: x0, y: y0 },
            lin_vel: Vec2 { x: vx0, y: vy0 },
            ang_pos: t0,
            ang_vel: o0,
            ..Default::default()
        }
    }
}
impl Default for Body {
    fn default() -> Self {
        Self {
            lin_pos: Vec2::new(0.0, 0.0),
            lin_vel: Vec2::new(0.0, 0.0),
            lin_acc: Vec2::new(0.0, 0.0),

            ang_pos: 0.0,
            ang_vel: 0.0,
            ang_acc: 0.0,
        }
    }
}

#[derive(Clone)]
struct Asteroid {
    body: Body,
    sides: u8,
    size: f32,
}

impl Asteroid {
    fn new(body: Body, sides: u8, size: f32) -> Self {
        Self { body, sides, size }
    }
}

trait Draw {
    fn draw(&self) -> ();
}

impl Draw for Asteroid {
    fn draw(&self) -> () {
        draw_poly_lines(
            self.body.lin_pos.x,
            self.body.lin_pos.y,
            self.sides,
            self.size,
            self.body.ang_pos.to_degrees(),
            6.0,
            WHITE,
        );
    }
}

trait Update {
    fn update(&mut self, dt: f32) -> ();
}

impl Update for Body {
    fn update(&mut self, dt: f32) -> () {
        self.lin_vel += self.lin_acc * dt;
        self.lin_pos += self.lin_vel * dt;

        self.ang_vel += self.ang_acc * dt;
        self.ang_pos += self.ang_vel * dt;
    }
}
