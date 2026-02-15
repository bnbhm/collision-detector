use macroquad::prelude::*;
use std::f32::consts::TAU;

#[macroquad::main("Braideroids : Asteroids = Braid;")]
async fn main() {
    let mut asteroid1 = Asteroid::new(Body::new(100.0, 100.0, 20.0, 50.0, 0.0, 2.0), 3, 120.0);
    let mut asteroid2 = Asteroid::new(Body::new(250.0, 100.0, 100.0, -20.0, 5.0, -0.4), 6, 50.0);

    let mut game_last_tick = get_time() as f32;
    loop {
        let current_tick = get_time() as f32;
        let dt = current_tick - game_last_tick;
        game_last_tick = current_tick;

        asteroid1.body.update(dt);
        asteroid2.body.update(dt);

        reflect(&mut asteroid1.body);
        reflect(&mut asteroid2.body);

        let collided = collision(&asteroid1, &asteroid2);

        asteroid1.draw();
        asteroid2.draw();
        if collided {
            draw_poly(
                asteroid1.body.lin_pos.x,
                asteroid1.body.lin_pos.y,
                asteroid1.sides,
                asteroid1.size,
                asteroid1.body.ang_pos.to_degrees(),
                WHITE,
            );
            draw_poly(
                asteroid2.body.lin_pos.x,
                asteroid2.body.lin_pos.y,
                asteroid2.sides,
                asteroid2.size,
                asteroid2.body.ang_pos.to_degrees(),
                WHITE,
            );
        }

        let vertices1 = asteroid1.shape();
        vertices1.iter().for_each(|vertice| {
            draw_circle(vertice.x, vertice.y, 5.0, BLUE);
        });

        let vertices2 = asteroid2.shape();
        vertices2.iter().for_each(|vertice| {
            draw_circle(vertice.x, vertice.y, 5.0, BLUE);
        });

        next_frame().await;
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

fn collision(object1: &impl Shape, object2: &impl Shape) -> bool {
    for vertice in object1.shape() {
        if inside(&vertice, object2) {
            return true;
        }
    }
    for vertice in object2.shape() {
        if inside(&vertice, object1) {
            return true;
        }
    }
    false
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
    if body.lin_pos.x >= screen_width() {
        body.lin_vel = Vec2::new(-1.0 * body.lin_vel.x, body.lin_vel.y)
    }
    if body.lin_pos.x <= 0.0 {
        body.lin_vel = Vec2::new(-1.0 * body.lin_vel.x, body.lin_vel.y)
    }
    if body.lin_pos.y >= screen_height() {
        body.lin_vel = Vec2::new(body.lin_vel.x, -1.0 * body.lin_vel.y)
    }
    if body.lin_pos.y <= 0.0 {
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
