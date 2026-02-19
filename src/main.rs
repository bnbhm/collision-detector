use macroquad::prelude::*;
use std::f32::consts::TAU;

#[macroquad::main("Braideroids : Asteroids = Braid;")]
async fn main() {
    let mut asteroid1 = Asteroid::new(Body::new(600.0, 500.0, -25.0, 0.0, -2.0, 0.5), 3, 120.0);
    let mut asteroid2 = Asteroid::new(Body::new(300.0, 500.0, 25.0, 0.0, 0.0, -1.0), 3, 50.0);
    let mut asteroid3 = Asteroid::new(Body::new(500.0, 600.0, 0.0, 25.0, 0.0, -0.3), 5, 50.0);
    let mut asteroid4 = Asteroid::new(Body::new(500.0, 300.0, 0.0, -25.0, 0.0, 0.5), 5, 100.0);

    let mut sentient_timer = get_time() as f32;
    let mut sentient_delta = 0.0;
    let mut vanish_timer = get_time() as f32;

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

        highlight_collision(&mut asteroid1, &mut asteroid2);
        highlight_collision(&mut asteroid1, &mut asteroid3);
        highlight_collision(&mut asteroid1, &mut asteroid4);

        highlight_collision(&mut asteroid2, &mut asteroid3);
        highlight_collision(&mut asteroid2, &mut asteroid4);

        highlight_collision(&mut asteroid3, &mut asteroid4);

        if (get_time() as f32 - sentient_timer) <= (5.0 + sentient_delta) {
            draw_text(
                "Sometimes the blocks behave sentient.",
                50.0,
                100.0,
                32.0,
                WHITE,
            );
            vanish_timer = get_time() as f32;
        } else {
            if (get_time() as f32 - vanish_timer) >= 5.0 {
                sentient_timer = get_time() as f32;
                vanish_timer = get_time() as f32;
                sentient_delta = 5.0 * (rand::rand() as f32 / u32::MAX as f32);
            }
        }

        next_frame().await;
    }
}

fn draw_vertices(object: &impl Shape) {
    let vertices1 = object.shape();
    vertices1.iter().for_each(|vertice| {
        draw_circle(vertice.x, vertice.y, 5.0, BLUE);
    });
}

fn highlight_collision(asteroid1: &mut Asteroid, asteroid2: &mut Asteroid) {
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
        collision_info[0].iter().for_each(|vertice_collision_info| {
            let vertice = asteroid1.shape()[vertice_collision_info.vertice_index];
            draw_circle(vertice.x, vertice.y, 4.0, RED);

            let vertice_p = asteroid1.shape()[(vertice_collision_info.vertice_index
                + asteroid1.shape().len()
                - 1)
                % asteroid1.shape().len()];
            let vertice_n = asteroid1.shape()
                [(vertice_collision_info.vertice_index + 1) % asteroid1.shape().len()];
            draw_line(vertice.x, vertice.y, vertice_p.x, vertice_p.y, 3.0, RED);
            draw_line(vertice.x, vertice.y, vertice_n.x, vertice_n.y, 3.0, RED);

            let mtv = 1.0
                * vertice_collision_info.mtv.min_overlap_value
                * vertice_collision_info.mtv.min_overlap_axis;

            draw_line(
                vertice.x,
                vertice.y,
                vertice.x + mtv.x,
                vertice.y + mtv.y,
                2.0,
                WHITE,
            );
            asteroid1.body.lin_acc = 100.0 * mtv;
            asteroid1.body.ang_acc = (vertice - asteroid1.body.lin_pos).perp_dot(mtv);
        });
        collision_info[1].iter().for_each(|vertice_collision_info| {
            let vertice = asteroid2.shape()[vertice_collision_info.vertice_index];
            draw_circle(vertice.x, vertice.y, 4.0, RED);

            let vertice_p = asteroid2.shape()[(vertice_collision_info.vertice_index
                + asteroid2.shape().len()
                - 1)
                % asteroid2.shape().len()];
            let vertice_n = asteroid2.shape()
                [(vertice_collision_info.vertice_index + 1) % asteroid2.shape().len()];
            draw_line(vertice.x, vertice.y, vertice_p.x, vertice_p.y, 3.0, RED);
            draw_line(vertice.x, vertice.y, vertice_n.x, vertice_n.y, 3.0, RED);

            let mtv = 1.0
                * vertice_collision_info.mtv.min_overlap_value
                * vertice_collision_info.mtv.min_overlap_axis;

            draw_line(
                vertice.x,
                vertice.y,
                vertice.x + mtv.x,
                vertice.y + mtv.y,
                2.0,
                WHITE,
            );

            asteroid2.body.lin_acc = 100.0 * mtv;
            asteroid2.body.ang_acc = (vertice - asteroid2.body.lin_pos).perp_dot(mtv);
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

struct VerticeCollisionInfo {
    vertice_index: usize,
    mtv: MTV,
}

fn collision(object1: &impl Shape, object2: &impl Shape) -> Option<[Vec<VerticeCollisionInfo>; 2]> {
    let mut info1 = Vec::new();
    let mut info2 = Vec::new();

    for (vertice_index, vertice) in object1.shape().iter().enumerate() {
        if let Some(mtv) = inside(&vertice, object2) {
            info1.push(VerticeCollisionInfo { vertice_index, mtv });
        }
    }

    for (vertice_index, vertice) in object2.shape().iter().enumerate() {
        if let Some(mtv) = inside(&vertice, object1) {
            info2.push(VerticeCollisionInfo { vertice_index, mtv });
        }
    }

    if info1.is_empty() && info2.is_empty() {
        return None;
    } else {
        return Some([info1, info2]);
    }
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

struct MTV {
    min_overlap_value: f32,
    min_overlap_axis: Vec2,
}

fn inside(vertice: &Vec2, object: &impl Shape) -> Option<MTV> {
    let vertices = object.shape();
    let mut surface_perps = vec![];

    for it in 0..vertices.len() {
        let b = vertices[(it + 1) % vertices.len()];
        let a = vertices[it % vertices.len()];
        let normal = b - a;
        let orthogonal = rotate(&normal, TAU / 4.0).normalize();

        surface_perps.push(orthogonal);
    }

    let mut min_overlap_value = f32::MAX;
    let mut min_overlap_axis = surface_perps[0];
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
            return None;
        }

        if vertice_proj > *max_proj {
            return None;
        }

        let overlap = (*max_proj - vertice_proj).min(vertice_proj - *min_proj);

        if overlap < min_overlap_value {
            min_overlap_value = overlap;
            min_overlap_axis = axis;
        }
    }

    return Some(MTV {
        min_overlap_value,
        min_overlap_axis,
    });
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
        if self.lin_vel.length() > 100.0 {
            self.lin_vel = 100.0 * self.lin_vel / self.lin_vel.length();
        }
        self.lin_pos += self.lin_vel * dt;

        self.ang_vel += self.ang_acc * dt;
        if self.ang_vel.abs() > 3.0 {
            self.ang_vel = 3.0 * self.ang_vel / self.ang_vel.abs();
        }
        self.ang_pos += self.ang_vel * dt;

        self.lin_acc = Vec2::new(0.0, 0.0);
        self.ang_acc = 0.0;
    }
}
