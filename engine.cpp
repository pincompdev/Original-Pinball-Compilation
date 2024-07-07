#include <SDL.h>
#include <SDL_image.h>
#include <SDL_mixer.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <cstdarg>

int DISPLAY_WIDTH = 640;
int DISPLAY_HEIGHT = 480;
int FRAMES_PER_SECOND = 60;
double MS_PER_FRAME = 1000.0 / FRAMES_PER_SECOND;
int PHYSICS_TICKS_PER_FRAME = 24 * 60 / FRAMES_PER_SECOND;
double SIM_TIME_PER_PHYSICS_TICK = 1.0 / (24.0 * 60.0); //1440 Hz physics regardless of framerate
bool FULLSCREEN = false;
bool BORDERLESS = false;
int FREQUENCY = 44100;
int OUTPUT_CHANNELS = 2;
int SAMPLE_SIZE = 4096;
int CHANNELS_PER_TRACK = 64;
int PIXELS_PER_INCH = 32;
double PI = 3.14159265;
bool DISABLE_FLASHING = false;
double ANALOG_INPUT_EPSILON = 0.05;
double UI_SCALE = 1.0;
double HUD_SCALE = 1.0;
double CAMERA_PAN_SENSITIVITY = 8.0;
double CAMERA_ZOOM_SENSITIVITY = 2.0;
double GLOBAL_FRICTION = 0.6; //[0.0, 1.0], affects solids only
std::string CONTROLS_PATH = "controls/controls.txt";

class Point
{
public:
    double x;
    double y;
    Point(double x, double y): x(x), y(y)
    {

    }
    double square_magnitude()
    {
        return x * x + y * y;
    }
    double magnitude()
    {
        return std::sqrt(x * x + y * y);
    }
    double dot_product(Point& p)
    {
        return x * p.x + y * p.y;
    }
};

class View
{
public:
    Point center = Point(10.125, 21.0);
    Point pending_center = Point(10.125, 21.0);
    double zoom = 1.0;
    double pending_zoom = 1.0;
    double transition_time = 0.0;
    int pixel_width = DISPLAY_WIDTH;
    int pixel_height = DISPLAY_HEIGHT;
    bool debug = true; //TODO make debug lines appear above sprites regardless of collision layers
    bool sprites = true;
    bool HUD = true;
    bool to_screen_space(Point* real, Point* screen)
    {
        //Returns true if point is on screen, false if not
        screen->x = (real->x - center.x) * zoom * PIXELS_PER_INCH + DISPLAY_WIDTH / 2;
        screen->y = (real->y - center.y) * zoom * PIXELS_PER_INCH + DISPLAY_HEIGHT / 2;
        return !(screen->x < 0 || screen->x >= pixel_width || screen->y < 0 || screen->y >= pixel_height);
    }
    bool bounding_box_visible(Point* top_left, Point* bottom_right)
    {
        Point top_left_pixel = Point(0.0, 0.0);
        Point bottom_right_pixel = Point(0.0, 0.0);
        to_screen_space(top_left, &top_left_pixel);
        to_screen_space(bottom_right, &bottom_right_pixel);
        return !(bottom_right_pixel.x < 0 || top_left_pixel.x >= DISPLAY_WIDTH || bottom_right_pixel.y < 0 || top_left_pixel.y >= DISPLAY_HEIGHT);
    }
    void pass_time(double ms)
    {
        if (transition_time == 0.0)
        {
            return;
        }
        if (ms > transition_time)
        {
            transition_time = 0.0;
            center.x = pending_center.x;
            center.y = pending_center.y;
            zoom = pending_zoom;
        } else
        {
            center.x = center.x + (pending_center.x - center.x) * ms / transition_time;
            center.y = center.y + (pending_center.y - center.y) * ms / transition_time;
            zoom = zoom + (pending_zoom - zoom) * ms / transition_time;
        }
        transition_time -= ms;
    }
};

class Sprite
{
public:
    SDL_Texture* texture;
    double origin_x = 0;
    double origin_y = 0;
    int width = 0;
    int height = 0;
    SDL_RendererFlip flip = SDL_FLIP_NONE;
    Sprite(SDL_Renderer* renderer, const char filename[], double origin_x, double origin_y, int width, int height) : origin_x(origin_x), origin_y(origin_y), width(width), height(height)
    {
        //Note: filename begins in .exe directory and must include any subdirectories explicitly
        //Note: transparency is inherited from the imported image
        SDL_Surface* surface = IMG_Load(filename);
        texture = SDL_CreateTextureFromSurface(renderer, surface);
        SDL_FreeSurface(surface);
    }
    ~Sprite()
    {
        SDL_DestroyTexture(texture);
    }
    bool is_culled(View* view, double x, double y, int display_w, int display_h)
    {
        Point real_point = Point(x, y);
        Point screen_point = Point(0.0, 0.0);
        view->to_screen_space(&real_point, &screen_point);
        double sprite_length = (width + height) * view->zoom;
        return screen_point.x + sprite_length < 0 || screen_point.y + sprite_length < 0 || screen_point.x - sprite_length > display_w || screen_point.y - sprite_length > display_h;
    }
    void render(SDL_Renderer* renderer, View* view, double x, double y, double angle, int animation, int frame)
    {
        if (!view->sprites || is_culled(view, x, y, DISPLAY_WIDTH, DISPLAY_HEIGHT))
        {
            return;
        }
        angle *= 180.0 / PI;
        SDL_Rect src;
        src.x = width * frame;
        src.y = 0;
        src.w = width;
        src.h = height;
        Point real_top_left = Point(x - origin_x / PIXELS_PER_INCH, y - origin_y / PIXELS_PER_INCH);
        Point screen_top_left = Point(0.0, 0.0);
        view->to_screen_space(&real_top_left, &screen_top_left);
        SDL_FRect dst;
        dst.x = screen_top_left.x;
        dst.y = screen_top_left.y;
        dst.w = width * view->zoom;
        dst.h = height * view->zoom;
        SDL_FPoint center;
        center.x = origin_x * view->zoom;
        center.y = origin_y * view->zoom;
        SDL_RenderCopyExF(renderer, texture, &src, &dst, angle, &center, flip);
    }
    void dmd_render(SDL_Renderer* renderer, double x, double y, double angle, int animation, int frame)
    {
        angle *= 180.0 / PI;
        SDL_Rect src;
        src.x = width * frame;
        src.y = height * animation;
        src.w = width;
        src.h = height;
        Point top_left = Point(x - origin_x, y - origin_y);
        SDL_FRect dst;
        dst.x = std::floor(top_left.x);
        dst.y = std::floor(top_left.y);
        dst.w = width;
        dst.h = height;
        SDL_FPoint center;
        center.x = origin_x;
        center.y = origin_y;
        SDL_RenderCopyExF(renderer, texture, &src, &dst, angle, &center, flip);
    }
};

enum AudioTrack
{
    TRACK_PHYSICAL,
    TRACK_TABLE,
    TRACK_MENU,
    TRACK_MUSIC,
    TRACK_COUNT
};

class Sound
{
public:
    Mix_Chunk* chunk = nullptr;
    AudioTrack track;
    int loop_channel = -1;
    Sound(const char filename[], AudioTrack track) : track(track)
    {
        chunk = Mix_LoadWAV(filename);
    }
    ~Sound()
    {
        if (loop_channel >= 0)
        {
            Mix_HaltChannel(loop_channel);
        }
        if (chunk)
        {
            Mix_FreeChunk(chunk);
        }
    }
    void play(int volume = 128)
    {
        int channel = Mix_GroupAvailable(track);
        if (channel < 0)
        {
            channel = Mix_AllocateChannels(1) - 1;
            Mix_GroupChannel(channel, track);
            channel = Mix_GroupAvailable(track);
            if (channel < 0)
            {
                channel = Mix_GroupOldest(track);
            }
        }
        Mix_Volume(channel, volume);
        Mix_PlayChannel(channel, chunk, 0);
    }
    void loop(int volume)
    {
        //Note: new Sound object should be instantiated for each concurrent looping sound, even if the source file is the same
        if (loop_channel < 0)
        {
            int channel = Mix_GroupAvailable(track);
            if (channel < 0)
            {
                channel = Mix_AllocateChannels(1) - 1;
                Mix_GroupChannel(channel, track);
                channel = Mix_GroupAvailable(track);
                if (channel < 0)
                {
                    channel = Mix_GroupOldest(track);
                }
            }
            loop_channel = channel;
        }
        Mix_Volume(loop_channel, volume);
        Mix_PlayChannel(loop_channel, chunk, -1);
    }
    void set_loop_volume(int volume)
    {
        if (loop_channel >= 0)
        {
            Mix_Volume(loop_channel, volume);
        }
    }
};

class Ball
{
public:
    bool tangible = true;
    Point position = Point(0.0, 0.0);
    Point velocity = Point(0.0, 0.0);
    double angular_displacement = 0.0;
    double angular_velocity = 0.0;
    double radius = 0.53125; //inches
    double mass = 2.6; //ounces
    double moment_ratio = 0.4;
    Sound* impact_sound = nullptr;
    double min_volume_energy = 432.166146;
    double max_volume_energy = 13829.3167;
    Sound* roll_sound = nullptr;
    double max_roll_volume_speed = 72.931255;
    double detection_speed = 0.0;
    double friction_multiplier = 1.0;
    int layer = 0;
    int layer_pending = 0; //used to send data from layer portal to script
    bool sprite_rotation = false;
    bool magnetic = true;
    bool show_spin = true;
    Ball(double x, double y)
    {
        position.x = x;
        position.y = y;
    }
    void add_roll_sound(Sound* sound)
    {
        roll_sound = sound;
        roll_sound->loop(128);
    }
    void bounce(Point* unit_normal, double bounce_coefficient, double friction_coefficient, Sound* impact_sound = nullptr, Point* detection_unit_normal = nullptr)
    {
        if (!unit_normal)
        {
            return;
        }
        if (velocity.dot_product(*unit_normal) >= 0.0)
        {
            return;
        }
        double impact_speed = -velocity.dot_product(*unit_normal);
        if (detection_unit_normal)
        {
            double normal_detection_alignment = unit_normal->dot_product(*detection_unit_normal);
            Point adjusted_normal = Point(unit_normal->x * normal_detection_alignment, unit_normal->y * normal_detection_alignment);
            detection_speed = -velocity.dot_product(adjusted_normal);
        } else
        {
            detection_speed = impact_speed;
        }
        if (impact_sound)
        {
            double energy = mass * impact_speed * impact_speed;
            if (energy >= min_volume_energy)
            {
                int volume = energy * 128 / max_volume_energy;
                if (volume > 128)
                {
                    volume = 128;
                }
                impact_sound->play(volume);
            }
        }
        Point orthogonal_unit_normal = Point(-unit_normal->y, unit_normal->x);
        double vx1 = velocity.dot_product(orthogonal_unit_normal);
        double A = 1.0 + 1.0 / moment_ratio;
        //Note: friction coefficient can be negative (-1.0 corresponds to frictionless surface)
        double vx2 = vx1 - (1.0 + friction_coefficient) * (vx1 - radius * angular_velocity) / A;
        double new_angular_velocity = angular_velocity + (vx1 - vx2) / (moment_ratio * radius);

        Point normal_velocity = Point(unit_normal->x * impact_speed * (bounce_coefficient - 1.0), unit_normal->y * impact_speed * (bounce_coefficient - 1.0));
        Point tangential_velocity = Point(orthogonal_unit_normal.x * vx2, orthogonal_unit_normal.y * vx2);
        Point outgoing_velocity = Point(normal_velocity.x + tangential_velocity.x, normal_velocity.y + tangential_velocity.y);

        velocity.x = outgoing_velocity.x;
        velocity.y = outgoing_velocity.y;
        angular_velocity = new_angular_velocity;
    }
    bool is_touching(Ball* b, bool clip, double bounce_factor)
    {
        if (this == b)
        {
            return false;
        }
        Point center_to_center = Point(b->position.x - position.x, b->position.y - position.y);
        double center_to_center_magnitude = center_to_center.magnitude();
        if (center_to_center_magnitude < radius + b->radius)
        {
            if (clip)
            {
                Point center_to_center_unit = Point(center_to_center.x / center_to_center_magnitude, center_to_center.y / center_to_center_magnitude);
                double clip_vector_magnitude = radius + b->radius - center_to_center_magnitude;
                Point clip_vector_a = Point(-center_to_center_unit.x * clip_vector_magnitude / 2, -center_to_center_unit.y * clip_vector_magnitude / 2);
                Point clip_vector_b = Point(center_to_center_unit.x * clip_vector_magnitude / 2, center_to_center_unit.y * clip_vector_magnitude / 2);
                position.x += clip_vector_a.x;
                position.y += clip_vector_a.y;
                b->position.x += clip_vector_b.x;
                b->position.y += clip_vector_b.y;
                if (bounce_factor)
                {
                    Point relative_velocity = Point(velocity.x - b->velocity.x, velocity.y - b->velocity.y);
                    double relative_velocity_magnitude = relative_velocity.dot_product(center_to_center_unit);
                    Point relative_normal_velocity = Point(center_to_center_unit.x * relative_velocity_magnitude, center_to_center_unit.y * relative_velocity_magnitude);
                    double impact_speed = relative_normal_velocity.magnitude() / 2.0;
                    if (impact_sound)
                    {
                        double energy = (mass + b->mass) * impact_speed * impact_speed / 2.0;
                        if (energy >= min_volume_energy)
                        {
                            int volume = energy * 128 / max_volume_energy;
                            if (volume > 128)
                            {
                                volume = 128;
                            }
                            impact_sound->play(volume);
                        }
                    }
                    if (b->impact_sound)
                    {
                        double energy = (mass + b->mass) * impact_speed * impact_speed / 2.0;
                        if (energy >= b->min_volume_energy)
                        {
                            int volume = energy * 128 / b->max_volume_energy;
                            if (volume > 128)
                            {
                                volume = 128;
                            }
                            b->impact_sound->play(volume);
                        }
                    }
                    velocity.x -= bounce_factor * relative_normal_velocity.x * b->mass / (mass + b->mass);
                    velocity.y -= bounce_factor * relative_normal_velocity.y * b->mass / (mass + b->mass);
                    b->velocity.x += bounce_factor * relative_normal_velocity.x * mass / (mass + b->mass);
                    b->velocity.y += bounce_factor * relative_normal_velocity.y * mass / (mass + b->mass);
                    double relative_angular_velocity = angular_velocity + b->angular_velocity;
                    angular_velocity -= relative_angular_velocity * (b->mass * b->moment_ratio) / (mass * moment_ratio + b->mass * b->moment_ratio);
                    b->angular_velocity -= relative_angular_velocity * (mass * moment_ratio) / (mass * moment_ratio + b->mass * b->moment_ratio);
                }
            }
            return true;
        }
        return false;
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        Point bounding_box_top_left = Point(position.x - radius, position.y - radius);
        Point bounding_box_bottom_right = Point(position.x + radius, position.y + radius);
        if (!view->bounding_box_visible(&bounding_box_top_left, &bounding_box_bottom_right))
        {
            return;
        }
        if (view->sprites && show_spin)
        {
            Point center = Point(PIXELS_PER_INCH * view->zoom * (position.x - view->center.x) + DISPLAY_WIDTH / 2, PIXELS_PER_INCH * view->zoom * (position.y - view->center.y) + DISPLAY_HEIGHT / 2);
            int marks = 3;
            for (int i = 0; i <= marks; ++i)
            {
                double angle = 2 * PI * i / marks + angular_displacement;
                Point offset = Point(std::cos(angle) * radius * PIXELS_PER_INCH * view->zoom, std::sin(angle) * radius * PIXELS_PER_INCH * view->zoom);
                SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
                double intensity = std::min(angular_velocity * angular_velocity / 1024.0, 1.0);
                SDL_SetRenderDrawColor(renderer, 255, 255, 255, 176 * intensity);
                SDL_RenderDrawLine(renderer, center.x + offset.x / 2.0, center.y + offset.y / 2.0, center.x + offset.x, center.y + offset.y);
                SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_NONE);
            }
        }
        if (view->debug)
        {
            //TODO signal physical properties on debug renderer (for all collider types)
            Point center = Point(PIXELS_PER_INCH * view->zoom * (position.x - view->center.x) + DISPLAY_WIDTH / 2, PIXELS_PER_INCH * view->zoom * (position.y - view->center.y) + DISPLAY_HEIGHT / 2);
            int precision = 32;
            Point last_offset = Point(0.0, 0.0);
            for (int i = 0; i <= precision; ++i)
            {
                double angle = 2 * PI * i / precision + angular_displacement;
                Point offset = Point(std::cos(angle) * radius * PIXELS_PER_INCH * view->zoom, std::sin(angle) * radius * PIXELS_PER_INCH * view->zoom);
                if (i)
                {
                    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
                    SDL_RenderDrawLine(renderer, center.x + last_offset.x, center.y + last_offset.y, center.x + offset.x, center.y + offset.y);
                } else
                {
                    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
                    SDL_RenderDrawLine(renderer, center.x, center.y, center.x + offset.x, center.y + offset.y);
                }
                last_offset = offset;
            }
        }
    }
};

void transform_point(Point* p, Point* origin, Point* offset, double angular_offset, Point* transformed)
{
    if (!(p && origin && offset && transformed))
    {
        *transformed = *p;
        return;
    }
    Point from_origin = Point(p->x - origin->x, p->y - origin->y);
    Point rotated_from_origin = Point(from_origin.x * cos(angular_offset) - from_origin.y * sin(angular_offset), from_origin.x * sin(angular_offset) + from_origin.y * cos(angular_offset));
    Point translated_from_origin = Point(rotated_from_origin.x + offset->x, rotated_from_origin.y + offset->y);
    transformed->x = origin->x + translated_from_origin.x;
    transformed->y = origin->y + translated_from_origin.y;
}

class Magnet
{
public:
    Point position = Point(0.0, 0.0);
    bool active = true;
    double magnetism = 0.0; //negative to repel
    double depth = 0.5;
    double radius = 0.0;
    bool hole = false;
    Magnet(double x, double y, double magnetism) : magnetism(magnetism)
    {
        position.x = x;
        position.y = y;
    }
    void attract(Ball* ball)
    {
        if (!active || (!ball->magnetic && !hole))
        {
            return;
        }
        Point from_ball = Point(position.x - ball->position.x, position.y - ball->position.y);
        if (from_ball.x > radius || from_ball.x < -radius || from_ball.y > radius || from_ball.y < -radius)
        {
            return;
        }
        double square_planar_distance = from_ball.square_magnitude();
        if (radius != 0.0 && square_planar_distance > radius * radius)
        {
            return;
        }
        double square_real_distance = square_planar_distance + depth * depth;
        double acceleration = magnetism / square_real_distance;
        double planar_distance = std::sqrt(square_planar_distance);
        if (planar_distance == 0.0)
        {
            return;
        }
        double angle = std::atan2(depth, planar_distance);
        Point from_ball_unit = Point(from_ball.x / planar_distance, from_ball.y / planar_distance);
        Point acceleration_vector = Point(from_ball_unit.x * acceleration * std::cos(angle), from_ball_unit.y * acceleration * std::cos(angle));
        ball->velocity.x += acceleration_vector.x * SIM_TIME_PER_PHYSICS_TICK;
        ball->velocity.y += acceleration_vector.y * SIM_TIME_PER_PHYSICS_TICK;
        if (!hole)
        {
            //TODO (Phase X) test behavior of hole magnets
            ball->friction_multiplier += acceleration * std::sin(angle) / 386.09;
        }
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        Point bounding_box_top_left = Point(position.x - radius, position.y - radius);
        Point bounding_box_bottom_right = Point(position.x + radius, position.y + radius);
        if (!view->bounding_box_visible(&bounding_box_top_left, &bounding_box_bottom_right))
        {
            return;
        }
        if (view->debug)
        {
            Point center = Point(PIXELS_PER_INCH * view->zoom * (position.x - view->center.x) + DISPLAY_WIDTH / 2, PIXELS_PER_INCH * view->zoom * (position.y - view->center.y) + DISPLAY_HEIGHT / 2);
            int precision = 32;
            Point last_offset = Point(0.0, 0.0);
            double pictured_radius = radius;
            if (!radius)
            {
                pictured_radius = 1.0;
            }
            for (int i = 0; i <= precision; ++i)
            {
                double angle = 2.0 * PI * i / precision;
                Point offset = Point(std::cos(angle) * pictured_radius * PIXELS_PER_INCH * view->zoom, std::sin(angle) * pictured_radius * PIXELS_PER_INCH * view->zoom);
                if (i)
                {
                    SDL_SetRenderDrawColor(renderer, 128, 128, 128, 255);
                    SDL_RenderDrawLine(renderer, center.x + last_offset.x, center.y + last_offset.y, center.x + offset.x, center.y + offset.y);
                }
                last_offset = offset;
            }
        }
    }
};

class Collider
{
public:
    Point bounding_box_top_left = Point(0.0, 0.0);
    Point bounding_box_bottom_right = Point(0.0, 0.0);
    virtual bool is_touching(Ball* b, bool clip, Point* unit_normal) = 0;
    virtual void render(SDL_Renderer* renderer, View* view) = 0;
    virtual void render(SDL_Renderer* renderer, View* view, Point* origin, Point* offset, double angular_offset) = 0;
    bool solid = true;
    bool mobile = false;
    int switch_id = -1;
    int gate = 0; //allows passage through one side (-1 = right-hand/outside, 0 = two-sided, 1 = left-hand/inside)
    //Note: a mobile collider can force a ball through a gate in the wrong direction
    //Note: layer portals should be gates and should only be accessible from the open side (include in documentation)
    bool portal = false; //signals script to send ball to destination_layer on collision
    int destination_layer = 0;
    double bounce_coefficient = 1.4;
    double friction_coefficient = 1.0;
    double sensitivity = 0.0;
    Sound* impact_sound = nullptr;
    bool has_detection_unit_normal = false;
    Point detection_unit_normal = Point(0.0, 0.0);
    int id = -1;
    virtual ~Collider()
    {

    }
};

class PointCollider : public Collider
{
public:
    Point position = Point(0.0, 0.0);
    PointCollider(double x, double y)
    {
        position.x = x;
        position.y = y;
        bounding_box_top_left.x = x;
        bounding_box_top_left.y = y;
        bounding_box_bottom_right.x = x;
        bounding_box_bottom_right.y = y;
    }
    bool is_touching(Ball* ball, bool clip, Point* unit_normal)
    {
        Point point_to_ball_center = Point(ball->position.x - position.x, ball->position.y - position.y);
        double square_magnitude = point_to_ball_center.square_magnitude();
        if (square_magnitude < ball->radius * ball->radius)
        {
            if (clip)
            {
                double magnitude = std::sqrt(square_magnitude);
                if (magnitude)
                {
                    double clip_magnitude = ball->radius - magnitude;
                    Point clip_vector = Point(point_to_ball_center.x / magnitude * clip_magnitude, point_to_ball_center.y / magnitude * clip_magnitude);
                    ball->position.x += clip_vector.x;
                    ball->position.y += clip_vector.y;
                }
                if (unit_normal)
                {
                    unit_normal->x = point_to_ball_center.x / magnitude;
                    unit_normal->y = point_to_ball_center.y / magnitude;
                }
            }
            return true;
        }
        return false;
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        render(renderer, view, nullptr, nullptr, 0.0);
    }
    void render(SDL_Renderer* renderer, View* view, Point* origin, Point* offset, double angular_offset)
    {
        if (!view->debug)
        {
            return;
        }
        if (!origin && !view->bounding_box_visible(&bounding_box_top_left, &bounding_box_bottom_right))
        {
            return;
        }
        Point transformed_point = Point(0.0, 0.0);
        transform_point(&position, origin, offset, angular_offset, &transformed_point);
        Point screen_point = Point(0.0, 0.0);
        view->to_screen_space(&transformed_point, &screen_point);
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderDrawPoint(renderer, screen_point.x, screen_point.y);
    }
};

class Spinner;

void rip_spinner(Spinner*, double);

class LineCollider : public Collider //Note: line colliders should have point colliders at both ends
{
public:
    Point a = Point(0.0, 0.0);
    Point b = Point(0.0, 0.0);
    Spinner* spinner = nullptr;
    bool band_bounce = false; //When true, collider is bouncier in the middle than at the ends
    bool outer_bounce_coefficient = 1.0;
    LineCollider(double x1, double y1, double x2, double y2)
    {
        a.x = x1;
        a.y = y1;
        b.x = x2;
        b.y = y2;
        bounding_box_top_left.x = x1 < x2 ? x1 : x2;
        bounding_box_top_left.y = y1 < y2 ? y1 : y2;
        bounding_box_bottom_right.x = x1 > x2 ? x1 : x2;
        bounding_box_bottom_right.y = y1 > y2 ? y1 : y2;
    }
    bool is_touching(Ball* ball, bool clip, Point* unit_normal)
    {
        Point a_to_b = Point(b.x - a.x, b.y - a.y);
        double a_to_b_magnitude = a_to_b.magnitude();
        Point a_to_b_unit = Point(a_to_b.x / a_to_b_magnitude, a_to_b.y / a_to_b_magnitude);
        Point a_to_ball = Point(ball->position.x - a.x, ball->position.y - a.y);
        double ball_projection = a_to_ball.dot_product(a_to_b_unit);
        if (ball_projection < 0.0 || ball_projection > a_to_b_magnitude)
        {
            return false;
        }
        Point a_to_projection = Point(ball_projection * a_to_b.x / a_to_b_magnitude, ball_projection * a_to_b.y / a_to_b_magnitude);
        Point line_to_ball = Point(a_to_ball.x - a_to_projection.x, a_to_ball.y - a_to_projection.y);
        double line_to_ball_square_magnitude = line_to_ball.square_magnitude();
        if (line_to_ball_square_magnitude < ball->radius * ball->radius)
        {
            if (spinner)
            {
                Point right_unit_normal = Point(-a_to_b_unit.y, a_to_b_unit.x);
                double contact_speed = -ball->velocity.dot_product(right_unit_normal);
                rip_spinner(spinner, contact_speed);
            }
            bool gate_pass = false;
            if (gate)
            {
                gate_pass = true;
                Point previous_ball_position = Point(ball->position.x - ball->velocity.x * SIM_TIME_PER_PHYSICS_TICK, ball->position.y - ball->velocity.y * SIM_TIME_PER_PHYSICS_TICK);
                Point a_to_previous_ball = Point(previous_ball_position.x - a.x, previous_ball_position.y - a.y);
                Point right_unit_normal = Point(-a_to_b_unit.y, a_to_b_unit.x);
                double distance_right_of_line = a_to_previous_ball.dot_product(right_unit_normal);
                double allowed_speed = ball->velocity.dot_product(right_unit_normal) * gate;
                double epsilon = ball->radius / 4.0;
                if (allowed_speed <= 0.0 && ((gate == -1 && distance_right_of_line <= -ball->radius + epsilon) || (gate == 1 && distance_right_of_line >= ball->radius - epsilon)))
                {
                    gate_pass = false;
                }
                if (gate_pass && distance_right_of_line * gate < 0.0)
                {
                    return false;
                }
            }
            if (clip && !gate_pass)
            {
                double line_to_ball_magnitude = std::sqrt(line_to_ball_square_magnitude);
                if (line_to_ball_magnitude != 0.0)
                {
                    Point line_to_ball_unit = Point(line_to_ball.x / line_to_ball_magnitude, line_to_ball.y / line_to_ball_magnitude);
                    Point clip_vector = Point((ball->radius - line_to_ball_magnitude) * line_to_ball_unit.x, (ball->radius - line_to_ball_magnitude) * line_to_ball_unit.y);
                    ball->position.x += clip_vector.x;
                    ball->position.y += clip_vector.y;
                    if (unit_normal)
                    {
                        unit_normal->x = line_to_ball_unit.x;
                        unit_normal->y = line_to_ball_unit.y;
                    }
                }
            }
            return true;
        }
        return false;
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        render(renderer, view, nullptr, nullptr, 0.0);
    }
    void render(SDL_Renderer* renderer, View* view, Point* origin, Point* offset, double angular_offset)
    {
        if (!view->debug)
        {
            return;
        }
        if (!origin && !view->bounding_box_visible(&bounding_box_top_left, &bounding_box_bottom_right))
        {
            return;
        }
        Point transformed_a = Point(0.0, 0.0);
        transform_point(&a, origin, offset, angular_offset, &transformed_a);
        Point transformed_b = Point(0.0, 0.0);
        transform_point(&b, origin, offset, angular_offset, &transformed_b);
        Point screen_a = Point(0.0, 0.0);
        view->to_screen_space(&transformed_a, &screen_a);
        Point screen_b = Point(0.0, 0.0);
        view->to_screen_space(&transformed_b, &screen_b);
        SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
        SDL_RenderDrawLine(renderer, screen_a.x, screen_a.y, screen_b.x, screen_b.y);
    }
};

double standardize_angle(double angle)
{
    return std::fmod(angle, 2 * PI);
}

class ArcCollider : public Collider //Note: arc colliders should have point colliders at both ends
{
public:
    Point center = Point(0.0, 0.0);
    double radius;
    double start_angle = 0.0;
    double end_angle = 2 * PI;
    ArcCollider(double x, double y, double radius) : radius(radius)
    {
        center.x = x;
        center.y = y;
        bounding_box_top_left.x = x - radius;
        bounding_box_top_left.y = y - radius;
        bounding_box_bottom_right.x = x + radius;
        bounding_box_bottom_right.y = y + radius;
    }
    ArcCollider(double x, double y, double radius, double start_angle, double end_angle) : radius(radius), start_angle(start_angle), end_angle(end_angle)
    {
        center.x = x;
        center.y = y;
        bounding_box_top_left.x = x - radius;
        bounding_box_top_left.y = y - radius;
        bounding_box_bottom_right.x = x + radius;
        bounding_box_bottom_right.y = y + radius;
        start_angle = standardize_angle(start_angle);
        end_angle = standardize_angle(end_angle);
        if (end_angle > start_angle)
        {
            std::swap(start_angle, end_angle);
        }
    }
    bool is_touching(Ball* ball, bool clip, Point* unit_normal)
    {
        Point center_to_ball = Point(ball->position.x - center.x, ball->position.y - center.y);
        if (end_angle - start_angle < 2 * PI)
        {
            double angle = std::atan2(center_to_ball.y, center_to_ball.x);
            if (!((angle - 2 * PI > start_angle && angle - 2 * PI < end_angle) || (angle > start_angle && angle < end_angle) || (angle + 2 * PI > start_angle && angle + 2 * PI < end_angle)))
            {
                return false;
            }
        }
        double center_to_ball_magnitude = center_to_ball.magnitude();
        if (center_to_ball_magnitude > radius && center_to_ball_magnitude - ball->radius < radius)
        {
            bool gate_pass = false;
            if (gate == -1)
            {
                gate_pass = true;
                Point center_to_ball_unit = Point(center_to_ball.x / center_to_ball_magnitude, center_to_ball.y / center_to_ball_magnitude);
                Point previous_ball_position = Point(ball->position.x - ball->velocity.x * SIM_TIME_PER_PHYSICS_TICK, ball->position.y - ball->velocity.y * SIM_TIME_PER_PHYSICS_TICK);
                Point center_to_previous_ball = Point(previous_ball_position.x - center.x, previous_ball_position.y - center.y);
                double distance_outside_arc = center_to_previous_ball.magnitude() - radius;
                double allowed_speed = -ball->velocity.dot_product(center_to_ball_unit);
                double epsilon = ball->radius / 4.0;
                if (allowed_speed <= 0.0 && distance_outside_arc <= -ball->radius + epsilon)
                {
                    gate_pass = false;
                }
                if (gate_pass)
                {
                    return false;
                }
            } else if (gate == 1)
            {
                gate_pass = true;
                Point center_to_ball_unit = Point(center_to_ball.x / center_to_ball_magnitude, center_to_ball.y / center_to_ball_magnitude);
                Point previous_ball_position = Point(ball->position.x - ball->velocity.x * SIM_TIME_PER_PHYSICS_TICK, ball->position.y - ball->velocity.y * SIM_TIME_PER_PHYSICS_TICK);
                Point center_to_previous_ball = Point(previous_ball_position.x - center.x, previous_ball_position.y - center.y);
                double distance_outside_arc = center_to_previous_ball.magnitude() - radius;
                double allowed_speed = ball->velocity.dot_product(center_to_ball_unit);
                double epsilon = ball->radius / 4.0;
                if (allowed_speed <= 0.0 &&  distance_outside_arc >= ball->radius - epsilon)
                {
                    gate_pass = false;
                }
            }
            if (clip && !gate_pass)
            {
                double clip_vector_magnitude = radius + ball->radius - center_to_ball_magnitude;
                Point clip_vector = Point(center_to_ball.x / center_to_ball_magnitude * clip_vector_magnitude, center_to_ball.y / center_to_ball_magnitude * clip_vector_magnitude);

                ball->position.x += clip_vector.x;
                ball->position.y += clip_vector.y;
                if (unit_normal)
                {
                    unit_normal->x = center_to_ball.x / center_to_ball_magnitude;
                    unit_normal->y = center_to_ball.y / center_to_ball_magnitude;
                }
            }
            return true;
        }
        if (center_to_ball_magnitude < radius && center_to_ball_magnitude + ball->radius > radius)
        {
            bool gate_pass = false;
            if (gate == -1)
            {
                gate_pass = true;
                Point center_to_ball_unit = Point(center_to_ball.x / center_to_ball_magnitude, center_to_ball.y / center_to_ball_magnitude);
                Point previous_ball_position = Point(ball->position.x - ball->velocity.x * SIM_TIME_PER_PHYSICS_TICK, ball->position.y - ball->velocity.y * SIM_TIME_PER_PHYSICS_TICK);
                Point center_to_previous_ball = Point(previous_ball_position.x - center.x, previous_ball_position.y - center.y);
                double distance_outside_arc = center_to_previous_ball.magnitude() - radius;
                double allowed_speed = -ball->velocity.dot_product(center_to_ball_unit);
                double epsilon = ball->radius / 4.0;
                if (allowed_speed <= 0.0 && distance_outside_arc <= -ball->radius + epsilon)
                {
                    gate_pass = false;
                }
            } else if (gate == 1)
            {
                gate_pass = true;
                Point center_to_ball_unit = Point(center_to_ball.x / center_to_ball_magnitude, center_to_ball.y / center_to_ball_magnitude);
                Point previous_ball_position = Point(ball->position.x - ball->velocity.x * SIM_TIME_PER_PHYSICS_TICK, ball->position.y - ball->velocity.y * SIM_TIME_PER_PHYSICS_TICK);
                Point center_to_previous_ball = Point(previous_ball_position.x - center.x, previous_ball_position.y - center.y);
                double distance_outside_arc = center_to_previous_ball.magnitude() - radius;
                double allowed_speed = ball->velocity.dot_product(center_to_ball_unit);
                double epsilon = ball->radius / 4.0;
                if (allowed_speed <= 0.0 &&  distance_outside_arc >= ball->radius - epsilon)
                {
                    gate_pass = false;
                }
                if (gate_pass)
                {
                    return false;
                }
            }
            double target_distance = radius - ball->radius;
            if (clip && center_to_ball_magnitude && target_distance && !gate_pass)
            {
                Point ball_to_center_unit = Point(-center_to_ball.x / center_to_ball_magnitude, -center_to_ball.y / center_to_ball_magnitude);
                double intersection_magnitude = ball->radius + center_to_ball_magnitude - radius;
                Point clip_vector = Point(ball_to_center_unit.x * intersection_magnitude, ball_to_center_unit.y * intersection_magnitude);

                ball->position.x += clip_vector.x;
                ball->position.y += clip_vector.y;
                if (unit_normal)
                {
                    unit_normal->x = ball_to_center_unit.x;
                    unit_normal->y = ball_to_center_unit.y;
                }
            }
            return true;
        }
        return false;
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        render(renderer, view, nullptr, nullptr, 0.0);
    }
    void render(SDL_Renderer* renderer, View* view, Point* mc_origin, Point* mc_offset, double mc_angular_offset)
    {
        if (!view->debug)
        {
            return;
        }
        if (!mc_origin && !view->bounding_box_visible(&bounding_box_top_left, &bounding_box_bottom_right))
        {
            return;
        }
        int precision = (end_angle - start_angle) / 0.196349 * radius * view->zoom * 3;
        Point last_offset = Point(0.0, 0.0);
        SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
        Point transformed_center = Point(0.0, 0.0);
        transform_point(&center, mc_origin, mc_offset, mc_angular_offset, &transformed_center);
        for (int i = 0; i <= precision; ++i)
        {
            double angle = start_angle + (end_angle - start_angle) * i / precision + mc_angular_offset;
            Point offset = Point(std::cos(angle) * radius, std::sin(angle) * radius);
            if (i)
            {
                SDL_RenderDrawLine(renderer, (transformed_center.x + last_offset.x - view->center.x) * PIXELS_PER_INCH * view->zoom + DISPLAY_WIDTH / 2, (transformed_center.y + last_offset.y - view->center.y) * PIXELS_PER_INCH * view->zoom + DISPLAY_HEIGHT / 2, (transformed_center.x + offset.x - view->center.x) * PIXELS_PER_INCH * view->zoom + DISPLAY_WIDTH / 2, (transformed_center.y + offset.y - view->center.y) * PIXELS_PER_INCH * view->zoom + DISPLAY_HEIGHT / 2);
            } else
            {
                SDL_RenderDrawPoint(renderer, (transformed_center.x + last_offset.x - view->center.x) * PIXELS_PER_INCH * view->zoom + DISPLAY_WIDTH / 2, (transformed_center.y + last_offset.y - view->center.y) * PIXELS_PER_INCH * view->zoom + DISPLAY_HEIGHT / 2);
            }
            last_offset = offset;
        }
    }
};

class CompositeCollider : public Collider
{
public:
    std::vector<ArcCollider*> arc_colliders;
    std::vector<LineCollider*> line_colliders;
    std::vector<PointCollider*> point_colliders;
    bool endpoints_added = false;
    bool is_touching(Ball* ball, bool clip, Point* unit_normal)
    {
        for (unsigned int i = 0; i < arc_colliders.size(); ++i)
        {
            if (arc_colliders[i]->is_touching(ball, clip, unit_normal))
            {
                return true;
            }
        }
        for (unsigned int i = 0; i < line_colliders.size(); ++i)
        {
            if (line_colliders[i]->is_touching(ball, clip, unit_normal))
            {
                return true;
            }
        }
        for (unsigned int i = 0; i < point_colliders.size(); ++i)
        {
            if (point_colliders[i]->is_touching(ball, clip, unit_normal))
            {
                return true;
            }
        }
        return false;
    }
    void update_bounding_box(Collider* c)
    {
        if (c->bounding_box_top_left.x < bounding_box_top_left.x)
        {
            bounding_box_top_left.x = c->bounding_box_top_left.x;
        }
        if (c->bounding_box_top_left.y < bounding_box_top_left.y)
        {
            bounding_box_top_left.y = c->bounding_box_top_left.y;
        }
        if (c->bounding_box_bottom_right.x > bounding_box_bottom_right.x)
        {
            bounding_box_bottom_right.x = c->bounding_box_bottom_right.x;
        }
        if (c->bounding_box_bottom_right.y > bounding_box_bottom_right.y)
        {
            bounding_box_bottom_right.y = c->bounding_box_bottom_right.y;
        }
    }
    void add_collider(ArcCollider* c)
    {
        update_bounding_box(c);
        arc_colliders.push_back(c);
    }
    void add_collider(LineCollider* c)
    {
        update_bounding_box(c);
        line_colliders.push_back(c);
    }
    void add_collider(PointCollider* c)
    {
        update_bounding_box(c);
        point_colliders.push_back(c);
    }
    void add_endpoints()
    {
        if (endpoints_added)
        {
            std::cout << "WARNING: redundant endpoints" << std::endl;
        }
        endpoints_added = true;
        for (unsigned int i = 0; i < arc_colliders.size(); ++i)
        {
            ArcCollider* c = arc_colliders[i];
            Point point_A = Point(c->center.x + c->radius * std::cos(c->start_angle), c->center.y + c->radius * std::sin(c->start_angle));
            Point point_B = Point(c->center.x + c->radius * std::cos(c->end_angle), c->center.y + c->radius * std::sin(c->end_angle));
            point_colliders.push_back(new PointCollider(point_A.x, point_A.y));
            point_colliders.push_back(new PointCollider(point_B.x, point_B.y));
        }
        for (unsigned int i = 0; i < line_colliders.size(); ++i)
        {
            LineCollider* c = line_colliders[i];
            Point point_A = Point(c->a.x, c->a.y);
            Point point_B = Point(c->b.x, c->b.y);
            point_colliders.push_back(new PointCollider(point_A.x, point_A.y));
            point_colliders.push_back(new PointCollider(point_B.x, point_B.y));
        }
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        render(renderer, view, nullptr, nullptr, 0.0);
    }
    void render(SDL_Renderer* renderer, View* view, Point* origin, Point* offset, double angular_offset)
    {
        if (!origin && !view->bounding_box_visible(&bounding_box_top_left, &bounding_box_bottom_right))
        {
            return;
        }
        for (unsigned int i = 0; i < arc_colliders.size(); ++i)
        {
            arc_colliders[i]->render(renderer, view, origin, offset, angular_offset);
        }
        for (unsigned int i = 0; i < line_colliders.size(); ++i)
        {
            line_colliders[i]->render(renderer, view, origin, offset, angular_offset);
        }
        for (unsigned int i = 0; i < point_colliders.size(); ++i)
        {
            point_colliders[i]->render(renderer, view, origin, offset, angular_offset);
        }
    }
};

class MobileCollider : public Collider
{
public:
    Collider* original_collider = nullptr;
    Point origin = Point(0.0, 0.0);
    Point offset = Point(0.0, 0.0);
    double angular_offset = 0.0;
    Point velocity = Point(0.0, 0.0);
    Point acceleration = Point(0.0, 0.0);
    double spring_constant = 0.0;
    double angular_velocity = 0.0;
    double angular_acceleration = 0.0;
    double angular_spring_constant = 0.0;
    bool bounce_enabled = true;
    Point translation_range_top_left = Point(0.0, 0.0);
    Point translation_range_bottom_right = Point(0.0, 0.0);
    double rotation_range_positive = 0.0;
    double rotation_range_negative = 0.0;
    double motion_bounce_coefficient = 1.0;
    double angular_motion_bounce_coefficient = 1.0;
    Sound* rotation_stop_sound = nullptr;
    double rotation_stop_min_volume_velocity = 0.05;
    double rotation_stop_max_volume_velocity = 4.0;
    Sound* translation_stop_sound = nullptr;
    double translation_stop_max_volume_velocity = 0.125;
    double translation_stop_min_volume_velocity = 4.0;
    std::vector<std::pair<Magnet*, Point>> synced_magnets;
    //Note: live catch elasticity should cause offset to cross back over flipper in approximately 1/3 the time it takes for the flipper to fire on a modern table
    bool live_catch_enabled = false;
    double live_catch_elasticity = 16384.0; //angular acceleration per radian of motion lag (equivalent to square Hertz)
    double live_catch_damping_coefficient = 0.5; //multiplies velocity down when offset crosses 0.0 (force travels through flipper)
    double live_catch_angular_offset = 0.0; //lags behind angular_offset according to elasticity
    double live_catch_angular_velocity = 0.0; //responsible for effect on ball
    MobileCollider(Collider* c, double origin_x, double origin_y, double left_range, double right_range, double up_range, double down_range, double cw_range, double ccw_range)
    {
        //Set both rotation ranges to 0.0 for unrestricted rotation
        mobile = true;
        original_collider = c;
        origin.x = origin_x;
        origin.y = origin_y;
        translation_range_top_left.x = origin.x - left_range;
        translation_range_top_left.y = origin.y - up_range;
        translation_range_bottom_right.x = origin.x + right_range;
        translation_range_bottom_right.y = origin.y + down_range;
        rotation_range_positive = cw_range;
        rotation_range_negative = -ccw_range;
        double radius = std::sqrt((c->bounding_box_bottom_right.x - c->bounding_box_top_left.x) * (c->bounding_box_bottom_right.x - c->bounding_box_top_left.x) + (c->bounding_box_bottom_right.y - c->bounding_box_top_left.y) * (c->bounding_box_bottom_right.y - c->bounding_box_top_left.y));
        bounding_box_top_left.x = c->bounding_box_top_left.x - left_range - radius;
        bounding_box_top_left.y = c->bounding_box_top_left.y - up_range - radius;
        bounding_box_bottom_right.x = c->bounding_box_bottom_right.x + right_range + radius;
        bounding_box_bottom_right.y = c->bounding_box_bottom_right.y + down_range + radius;
    }
    void sync_magnet(Magnet* magnet, double x, double y)
    {
        synced_magnets.push_back(std::pair<Magnet*, Point>(magnet, Point(origin.x + x, origin.y + y)));
    }
    void simulate()
    {
        velocity.x += (acceleration.x + spring_constant * -offset.x) * SIM_TIME_PER_PHYSICS_TICK;
        velocity.y += (acceleration.y + spring_constant * -offset.y) * SIM_TIME_PER_PHYSICS_TICK;
        offset.x += velocity.x * SIM_TIME_PER_PHYSICS_TICK;
        offset.y += velocity.y * SIM_TIME_PER_PHYSICS_TICK;
        if (origin.x + offset.x < translation_range_top_left.x)
        {
            if (translation_stop_sound)
            {
                double energy = velocity.x * velocity.x;
                if (energy >= translation_stop_min_volume_velocity * translation_stop_min_volume_velocity)
                {
                    int volume = energy * 128 / translation_stop_max_volume_velocity / translation_stop_max_volume_velocity;
                    if (volume > 128)
                    {
                        volume = 128;
                    }
                    impact_sound->play(volume);
                }
            }
            offset.x = translation_range_top_left.x - origin.x;
            if (bounce_enabled)
            {
                velocity.x -= velocity.x * motion_bounce_coefficient;
            } else
            {
                velocity.x = 0.0;
            }
        }
        if (origin.x + offset.x > translation_range_bottom_right.x)
        {
            if (translation_stop_sound)
            {
                double energy = velocity.x * velocity.x;
                if (energy >= translation_stop_min_volume_velocity * translation_stop_min_volume_velocity)
                {
                    int volume = energy * 128 / translation_stop_max_volume_velocity / translation_stop_max_volume_velocity;
                    if (volume > 128)
                    {
                        volume = 128;
                    }
                    impact_sound->play(volume);
                }
            }
            offset.x = translation_range_bottom_right.x - origin.x;
            if (bounce_enabled)
            {
                velocity.x -= velocity.x * motion_bounce_coefficient;
            } else
            {
                velocity.x = 0.0;
            }
        }
        if (origin.y + offset.y < translation_range_top_left.y)
        {
            if (translation_stop_sound)
            {
                double energy = velocity.y * velocity.y;
                if (energy >= translation_stop_min_volume_velocity * translation_stop_min_volume_velocity)
                {
                    int volume = energy * 128 / translation_stop_max_volume_velocity / translation_stop_max_volume_velocity;
                    if (volume > 128)
                    {
                        volume = 128;
                    }
                    impact_sound->play(volume);
                }
            }
            offset.y = translation_range_top_left.y - origin.y;
            if (bounce_enabled)
            {
                velocity.y -= velocity.y * motion_bounce_coefficient;
            } else
            {
                velocity.y = 0.0;
            }
        }
        if (origin.y + offset.y > translation_range_bottom_right.y)
        {
            if (translation_stop_sound)
            {
                double energy = velocity.y * velocity.y;
                if (energy >= translation_stop_min_volume_velocity * translation_stop_min_volume_velocity)
                {
                    int volume = energy * 128 / translation_stop_max_volume_velocity / translation_stop_max_volume_velocity;
                    if (volume > 128)
                    {
                        volume = 128;
                    }
                    impact_sound->play(volume);
                }
            }
            offset.y = translation_range_bottom_right.y - origin.y;
            if (bounce_enabled)
            {
                velocity.y -= velocity.y * motion_bounce_coefficient;
            } else
            {
                velocity.y = 0.0;
            }
        }
        angular_velocity += (angular_acceleration + angular_spring_constant * -angular_offset) * SIM_TIME_PER_PHYSICS_TICK;
        angular_offset += angular_velocity * SIM_TIME_PER_PHYSICS_TICK;
        if (rotation_range_positive || rotation_range_negative)
        {
            if (angular_offset > rotation_range_positive)
            {
                angular_offset = rotation_range_positive;
                if (bounce_enabled)
                {
                    angular_velocity -= angular_velocity * angular_motion_bounce_coefficient;
                } else
                {
                    angular_velocity = 0.0;
                }
                if (rotation_stop_sound)
                {
                    double energy = angular_velocity * angular_velocity;
                    if (energy >= rotation_stop_min_volume_velocity * rotation_stop_min_volume_velocity)
                    {
                        int volume = energy * 128 / rotation_stop_max_volume_velocity / rotation_stop_max_volume_velocity;
                        if (volume > 128)
                        {
                            volume = 128;
                        }
                        impact_sound->play(volume); //TODO play sound for translational motion bounces?
                    }
                }
            }
            if (angular_offset < rotation_range_negative)
            {
                angular_offset = rotation_range_negative;
                if (bounce_enabled)
                {
                    angular_velocity -= angular_velocity * angular_motion_bounce_coefficient;
                } else
                {
                    angular_velocity = 0.0;
                }
                if (rotation_stop_sound)
                {
                    double energy = angular_velocity * angular_velocity;
                    if (energy >= rotation_stop_min_volume_velocity * rotation_stop_min_volume_velocity)
                    {
                        int volume = energy * 128 / rotation_stop_max_volume_velocity / rotation_stop_max_volume_velocity;
                        if (volume > 128)
                        {
                            volume = 128;
                        }
                        impact_sound->play(volume);
                    }
                }
            }
        }
        for (unsigned int i = 0; i < synced_magnets.size(); ++i)
        {
            Magnet* magnet = synced_magnets[i].first;
            Point relative = synced_magnets[i].second;
            Point real = Point(0.0, 0.0);
            relative_to_real(&relative, &real);
            magnet->position.x = real.x;
            magnet->position.y = real.y;
        }
        live_catch_angular_velocity -= live_catch_elasticity * (live_catch_angular_offset - angular_offset) * SIM_TIME_PER_PHYSICS_TICK;
        live_catch_angular_offset += live_catch_angular_velocity * SIM_TIME_PER_PHYSICS_TICK;
        if ((live_catch_angular_offset - angular_offset) * live_catch_angular_velocity * SIM_TIME_PER_PHYSICS_TICK > 0.0 && std::abs(live_catch_angular_offset - angular_offset) < std::abs(live_catch_angular_velocity * SIM_TIME_PER_PHYSICS_TICK))
        {
            live_catch_angular_velocity *= live_catch_damping_coefficient;
        }
    }
    void real_to_relative(Point* real, Point* relative)
    {
        Point real_origin = Point(origin.x + offset.x, origin.y + offset.y);
        Point from_real_origin = Point(real->x - real_origin.x, real->y - real_origin.y);
        Point from_relative_origin = Point(from_real_origin.x * cos(-angular_offset) - from_real_origin.y * sin(-angular_offset), from_real_origin.x * sin(-angular_offset) + from_real_origin.y * cos(-angular_offset));
        relative->x = origin.x + from_relative_origin.x;
        relative->y = origin.y + from_relative_origin.y;
    }
    void relative_to_real(Point* relative, Point* real)
    {
        Point from_origin = Point(relative->x - origin.x, relative->y - origin.y);
        Point rotated_from_origin = Point(from_origin.x * cos(angular_offset) - from_origin.y * sin(angular_offset), from_origin.x * sin(angular_offset) + from_origin.y * cos(angular_offset));
        real->x = origin.x + rotated_from_origin.x + offset.x;
        real->y = origin.y + rotated_from_origin.y + offset.y;
    }
    void velocity_at_point(Point* p, Point* v)
    {
        Point from_origin = Point(p->x - origin.x - offset.x, p->y - origin.y - offset.y);
        double from_origin_magnitude = from_origin.magnitude();
        if (from_origin_magnitude == 0.0)
        {
            v->x = 0.0;
            v->y = 0.0;
            return;
        }
        Point from_origin_unit = Point(from_origin.x / from_origin_magnitude, from_origin.y / from_origin_magnitude);
        Point orthogonal_from_origin_unit = Point(-from_origin_unit.y, from_origin_unit.x);
        double speed_from_rotation = from_origin_magnitude * angular_velocity;
        Point velocity_from_rotation = Point(orthogonal_from_origin_unit.x * speed_from_rotation, orthogonal_from_origin_unit.y * speed_from_rotation);
        v->x = velocity_from_rotation.x + velocity.x;
        v->y = velocity_from_rotation.y + velocity.y;
    }
    void velocity_at_relative_point(Point* p, Point* v)
    {
        Point from_origin = Point(p->x - origin.x, p->y - origin.y);
        double from_origin_magnitude = from_origin.magnitude();
        if (from_origin_magnitude == 0.0)
        {
            v->x = velocity.x;
            v->y = velocity.y;
            return;
        }
        Point from_origin_unit = Point(from_origin.x / from_origin_magnitude, from_origin.y / from_origin_magnitude);
        Point orthogonal_from_origin_unit = Point(-from_origin_unit.y, from_origin_unit.x);
        double speed_from_rotation = from_origin_magnitude * angular_velocity;
        Point velocity_from_rotation = Point(orthogonal_from_origin_unit.x * speed_from_rotation, orthogonal_from_origin_unit.y * speed_from_rotation);
        Point rotated_velocity = Point(velocity.x * std::cos(-angular_offset) - velocity.y * std::sin(-angular_offset), velocity.x * std::sin(-angular_offset) + velocity.y * std::cos(angular_offset));
        v->x = velocity_from_rotation.x + rotated_velocity.x;
        v->y = velocity_from_rotation.y + rotated_velocity.y;
    }
    bool is_touching(Ball* b, bool clip, Point* unit_normal)
    {
        Point relative_ball_position = Point(0.0, 0.0);
        real_to_relative(&(b->position), &relative_ball_position);
        Ball relative_ball(relative_ball_position.x, relative_ball_position.y);
        relative_ball.radius = b->radius;
        relative_ball.tangible = b->tangible;
        Point relative_unit_normal = Point(0.0, 0.0);
        Point reference_velocity = Point(0.0, 0.0);
        Point relative_ball_position_after_clip = Point(relative_ball.position.x, relative_ball.position.y);
        Point relative_contact_point = Point(relative_ball_position_after_clip.x , relative_ball_position_after_clip.y);
        velocity_at_relative_point(&relative_contact_point, &reference_velocity);
        relative_ball.velocity.x = b->velocity.x * std::cos(-angular_offset) - b->velocity.y * std::sin(-angular_offset) - reference_velocity.x;
        relative_ball.velocity.y = b->velocity.x * std::sin(-angular_offset) + b->velocity.y * std::cos(-angular_offset) - reference_velocity.y;
        relative_ball.angular_velocity = b->angular_velocity;
        bool touching = original_collider->is_touching(&relative_ball, clip, &relative_unit_normal);
        if (!touching)
        {
            return false;
        }
        if (clip)
        {
            Point reference_velocity = Point(0.0, 0.0);
            Point relative_ball_position_after_clip = Point(relative_ball.position.x, relative_ball.position.y);
            Point relative_contact_point = Point(relative_ball_position_after_clip.x - relative_unit_normal.x * relative_ball.radius, relative_ball_position_after_clip.y - relative_unit_normal.y * relative_ball.radius);
            velocity_at_relative_point(&relative_contact_point, &reference_velocity);
            //1. adjust relative ball velocity according to movement of frame of reference
            relative_ball.velocity.x = b->velocity.x * std::cos(-angular_offset) - b->velocity.y * std::sin(-angular_offset) - reference_velocity.x;
            relative_ball.velocity.y = b->velocity.x * std::sin(-angular_offset) + b->velocity.y * std::cos(-angular_offset) - reference_velocity.y;
            double relative_live_catch_angular_offset = live_catch_angular_offset - angular_offset;
            double effective_bounce_coefficient = bounce_coefficient;
            double effective_friction_coefficient = friction_coefficient;
            Point orthogonal_origin_to_relative_contact_point = Point(origin.y - relative_contact_point.y, relative_contact_point.x - origin.x);
            if (live_catch_enabled && orthogonal_origin_to_relative_contact_point.dot_product(relative_unit_normal) * relative_live_catch_angular_offset > 0.0 && orthogonal_origin_to_relative_contact_point.dot_product(relative_unit_normal) * live_catch_angular_velocity > 0.0)
            {
                double live_catch_effectiveness = 1.0 - std::pow(1.35914091, -orthogonal_origin_to_relative_contact_point.magnitude() * std::abs(live_catch_angular_velocity)); //Note: e / 2 is selected arbitrarily as a reasonable base
                effective_bounce_coefficient = (bounce_coefficient - 1.0) * (1.0 - live_catch_effectiveness) + 1.0;
                //Code for live catch to increase friction
                //effective_friction_coefficient = (1.0 - (1.0 - (friction_coefficient + 1.0) / 2.0) * (1.0 - live_catch_effectiveness)) * 2.0 - 1.0;
                //Code for live catch to decrease friction
                effective_friction_coefficient = (friction_coefficient + 1.0) * (1.0 - live_catch_effectiveness) - 1.0;
            }
            //2. obtain new relative ball velocity after bouncing
            //TODO add relative detection unit normal
            relative_ball.bounce(&relative_unit_normal, effective_bounce_coefficient, effective_friction_coefficient, impact_sound);
            Point relative_clip_vector = Point(relative_ball_position_after_clip.x - relative_ball_position.x, relative_ball_position_after_clip.y - relative_ball_position.y);
            Point clip_vector = Point(relative_clip_vector.x * std::cos(angular_offset) - relative_clip_vector.y * std::sin(angular_offset), relative_clip_vector.x * std::sin(angular_offset) + relative_clip_vector.y * std::cos(angular_offset));
            b->position.x += clip_vector.x;
            b->position.y += clip_vector.y;
            //3. set real ball velocity by converting relative ball velocity
            b->velocity.x = relative_ball.velocity.x * std::cos(angular_offset) - relative_ball.velocity.y * std::sin(angular_offset);
            b->velocity.y = relative_ball.velocity.x * std::sin(angular_offset) + relative_ball.velocity.y * std::cos(angular_offset);
            b->angular_velocity = relative_ball.angular_velocity;
            //4. add difference in frame of reference
            Point real_reference_velocity = Point(0.0, 0.0);
            Point real_unit_normal = Point(relative_unit_normal.x * std::cos(angular_offset) - relative_unit_normal.y * std::sin(angular_offset), relative_unit_normal.x * std::sin(angular_offset) + relative_unit_normal.y * std::cos(angular_offset));
            Point real_contact_point = Point(b->position.x - real_unit_normal.x * b->radius, b->position.y - real_unit_normal.y * b->radius);
            velocity_at_point(&real_contact_point, &real_reference_velocity);
            b->velocity.x += real_reference_velocity.x;
            b->velocity.y += real_reference_velocity.y;
        }
        if (unit_normal)
        {
            unit_normal->x = relative_unit_normal.x * std::cos(angular_offset) - relative_unit_normal.y * std::sin(angular_offset);
            unit_normal->y = relative_unit_normal.x * std::sin(angular_offset) + relative_unit_normal.y * std::cos(angular_offset);
        }
        return true;
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        if (!view->bounding_box_visible(&bounding_box_top_left, &bounding_box_bottom_right))
        {
            return;
        }
        original_collider->render(renderer, view, &origin, &offset, angular_offset);
        if (live_catch_enabled)
        {
            original_collider->render(renderer, view, &origin, &offset, live_catch_angular_offset);
        }
    }
    void render(SDL_Renderer* renderer, View* view, Point* origin, Point* offset, double angular_offset)
    {
        //nested mobile colliders not implemented
    }
};

class MobileCollider;

class SpriteInstance
{
public:
    Sprite* sprite = nullptr;
    Point position = Point(0.0, 0.0);
    double angle = 0.0;
    int animation = 0;
    int frame = 0;
    double ms_elapsed_in_animation = 0;
    double animation_fps = 24.0;
    bool can_rotate = true;
    bool visible = true;
    int loop_frame_count = 0; //0 for static image
    int loops = -1; //-1 to loop infinitely
    int stop_frame = 0;
    bool paused = false;
    int loop_start_frame = 0;
    MobileCollider* synced_collider = nullptr;
    SpriteInstance(Sprite* sprite, double x = 0.0, double y = 0.0, double angle = 0.0) : sprite(sprite), angle(angle)
    {
        position.x = x;
        position.y = y;
    }
    void sync(Ball* b)
    {
        position.x = b->position.x;
        position.y = b->position.y;
        if (can_rotate && b->sprite_rotation)
        {
            angle = b->angular_displacement;
        }
    }
    void sync(MobileCollider* mc)
    {
        position.x = mc->origin.x + mc->offset.x;
        position.y = mc->origin.y + mc->offset.y;
        if (can_rotate)
        {
            angle = mc->angular_offset;
        }
    }
    void pass_time(double ms)
    {
        if (!loop_frame_count)
        {
            return;
        }
        if (!paused)
        {
            ms_elapsed_in_animation += ms;
        }
        double ms_per_frame = 1000.0 / animation_fps;
        frame = static_cast<int>(ms_elapsed_in_animation / ms_per_frame) % loop_frame_count + loop_start_frame;
        if (loops != -1 && ms_elapsed_in_animation / ms_per_frame / loop_frame_count > loops)
        {
            stop_animation();
        }
    }
    void start_animation(int animation_id, int frame_count, double frames_per_second, int loop_count = -1, int frame_to_stop = 0, int offset = 0, int start_frame = 0)
    {
        paused = false;
        ms_elapsed_in_animation = 0;
        animation = animation_id;
        frame = offset;
        loop_frame_count = frame_count;
        animation_fps = frames_per_second;
        loops = loop_count;
        stop_frame = frame_to_stop;
        loop_start_frame = start_frame;
    }
    void start_synchronous_animation(int animation_id, int frame_count, double frames_per_second, double sync_ms, int start_frame)
    {
        paused = false;
        ms_elapsed_in_animation = sync_ms;
        animation = animation_id;
        frame = std::floor(frames_per_second / sync_ms);
        loop_frame_count = frame_count;
        animation_fps = frames_per_second;
        loops = -1;
        stop_frame = 0;
        loop_start_frame = start_frame;
    }
    void stop_animation()
    {
        ms_elapsed_in_animation = 0;
        loop_frame_count = 0;
        frame = stop_frame;
        loop_start_frame = 0;
    }
    void pause_animation()
    {
        paused = true;
    }
    void resume_animation()
    {
        paused = false;
    }
    void toggle_animation()
    {
        paused = !paused;
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        if (!visible)
        {
            return;
        }
        sprite->render(renderer, view, position.x, position.y, angle, animation, frame);
    }
    void dmd_render(SDL_Renderer* renderer)
    {
        if (!visible)
        {
            return;
        }
        sprite->dmd_render(renderer, position.x, position.y, angle, animation, frame);
    }
};

enum TextAlignment
{
    ALIGN_BOTTOM_LEFT,
    ALIGN_BOTTOM,
    ALIGN_BOTTOM_RIGHT,
    ALIGN_LEFT,
    ALIGN_CENTER,
    ALIGN_RIGHT,
    ALIGN_TOP_LEFT,
    ALIGN_TOP,
    ALIGN_TOP_RIGHT
};

class FontSize
{
public:
    int height;
    int max_width;
    int width_per_char[128];
    int kerning_space = 1;
    SDL_Texture* char_set = nullptr;
    FontSize(SDL_Renderer* renderer, const char filename[], const char kerning_data_filename[])
    {
        SDL_Surface* surface = IMG_Load(filename);
        char_set = SDL_CreateTextureFromSurface(renderer, surface);
        SDL_FreeSurface(surface);
        std::ifstream kerning_data(kerning_data_filename);
        int i = 0;
        std::string line;
        while (i < 128)
        {
            std::getline(kerning_data, line);
            width_per_char[i] = std::stoi(line);
            ++i;
        }
        std::getline(kerning_data, line);
        kerning_space = std::stoi(line);
        std::getline(kerning_data, line);
        height = std::stoi(line);
        std::getline(kerning_data, line);
        max_width = std::stoi(line);
    }
    FontSize(SDL_Renderer* renderer, const char filename[], int monospaced_width, int height) : height(height)
    {
        SDL_Surface* surface = IMG_Load(filename);
        char_set = SDL_CreateTextureFromSurface(renderer, surface);
        SDL_FreeSurface(surface);
        max_width = monospaced_width;
        for (int i = 0; i < 128; ++i)
        {
            width_per_char[i] = monospaced_width;
        }
    }
    ~FontSize()
    {
        SDL_DestroyTexture(char_set);
    }
    int get_width(std::string text)
    {
        int total_width = 0;
        for (unsigned int i = 0; i < text.size(); ++i)
        {
            int c = text[i];
            total_width += width_per_char[c] + kerning_space;
        }
        return total_width;
    }
    void show_text(SDL_Renderer* renderer, std::string text, int x, int y, TextAlignment align = ALIGN_TOP_LEFT, double scale = 1.0)
    {
        switch (align)
        {
        case ALIGN_BOTTOM_LEFT:
            y -= height;
            break;
        case ALIGN_BOTTOM:
            x -= get_width(text) / 2;
            y -= height;
            break;
        case ALIGN_BOTTOM_RIGHT:
            x -= get_width(text);
            y -= height;
            break;
        case ALIGN_LEFT:
            y -= height / 2;
            break;
        case ALIGN_CENTER:
            x -= get_width(text) / 2;
            y -= height / 2;
            break;
        case ALIGN_RIGHT:
            x -= get_width(text);
            y -= height / 2;
            break;
        case ALIGN_TOP_LEFT:
            //Already aligned
            break;
        case ALIGN_TOP:
            x -= get_width(text) / 2;
            break;
        case ALIGN_TOP_RIGHT:
            x -= get_width(text);
            break;
        }
        SDL_FRect dst;
        dst.x = x;
        dst.y = y;
        dst.h = height * scale;
        SDL_Rect src;
        src.y = 0;
        src.h = height;
        for (unsigned int i = 0; i < text.size(); ++i)
        {
            int c = text[i];
            src.x = max_width * c;
            src.w = max_width;
            dst.w = max_width * scale;
            SDL_RenderCopyF(renderer, char_set, &src, &dst);
            dst.x += (width_per_char[c] + kerning_space) * scale;
        }
    }
};

class Font
{
public:
    std::deque<FontSize*> font_sizes;
    Font()
    {

    }
    void add_size(FontSize* font_size)
    {
        /*
        font_sizes.push_back(font_size);
        std::sort(font_sizes.begin(), font_sizes.end(), [](FontSize* a, FontSize* b)
        {
            if (a->height > b->height)
            {
                return true;
            } else if (a->height == b->height)
            {
                return a->max_width >= b->max_width;
            } else
            {
                return false;
            }
        });
        */
        //TODO (Phase X) test binary search (a and b may need to be swapped)
        int high = font_sizes.size();
        int low = 0;
        while (high > low)
        {
            int middle = (high + low) / 2;
            FontSize* a = font_sizes[middle];
            FontSize* b = font_size;
            if (a->height > b->height || (a->height == b->height && a->max_width >= b->max_width))
            {
                low = middle;
            } else
            {
                high = middle;
            }
        }
        font_sizes.insert(font_sizes.begin() + low, font_size);
    }
    void show_text(SDL_Renderer* renderer, std::string text, int x, int y, TextAlignment align = ALIGN_TOP_LEFT, int max_width = 0, int max_height = 0, double scale = 1.0)
    {
        //TODO (Phase X) test with multiple sizes
        FontSize* fit = nullptr;
        for (unsigned int i = 0; i < font_sizes.size(); ++i)
        {
            fit = font_sizes[i];
            if (max_height && fit->height * scale > max_height)
            {
                continue;
            }
            if (max_width && fit->get_width(text) * scale > max_width)
            {
                continue;
            }
            break;
        }
        fit->show_text(renderer, text, x, y, align, scale);
    }
};

class DotMatrixDisplay
{
public:
    int width = 96;
    int height = 32;
    Sprite* dot = nullptr;
    int unlit = 16;
    int lit_r = 255;
    int lit_g = 128;
    int lit_b = 0;
    bool invert = false;
    SDL_Texture* dot_matrix = nullptr;
    SDL_Texture* intermediate = nullptr;
    SDL_Texture* alpha_mask = nullptr;
    std::vector<SpriteInstance*> sprite_instances;
    DotMatrixDisplay(SDL_Renderer* renderer, int width, int height, Sprite* dot) : width(width), height(height), dot(dot)
    {
        dot_matrix = SDL_CreateTexture(renderer, 0, SDL_TEXTUREACCESS_TARGET, width, height);
        SDL_SetTextureBlendMode(dot_matrix, SDL_BLENDMODE_BLEND);
        SDL_SetRenderTarget(renderer, dot_matrix);
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        SDL_SetRenderTarget(renderer, nullptr);

        intermediate = SDL_CreateTexture(renderer, 0, SDL_TEXTUREACCESS_TARGET, width * dot->width, height * dot->height);
        SDL_SetTextureBlendMode(intermediate, SDL_BLENDMODE_BLEND);

        alpha_mask = SDL_CreateTexture(renderer, 0, SDL_TEXTUREACCESS_TARGET, width * dot->width, height * dot->height);
        generate_alpha_mask(renderer);
        SDL_SetTextureBlendMode(alpha_mask, SDL_ComposeCustomBlendMode(SDL_BLENDFACTOR_ZERO, SDL_BLENDFACTOR_ONE, SDL_BLENDOPERATION_ADD, SDL_BLENDFACTOR_ONE, SDL_BLENDFACTOR_ZERO, SDL_BLENDOPERATION_ADD));
        SDL_SetRenderTarget(renderer, nullptr);
    }
    ~DotMatrixDisplay()
    {
        SDL_DestroyTexture(dot_matrix);
    }
    void clear_display(SDL_Renderer* renderer, int brightness = 0)
    {
        SDL_SetRenderTarget(renderer, dot_matrix);
        SDL_SetRenderDrawColor(renderer, brightness, brightness, brightness, 255);
        SDL_RenderClear(renderer);
        SDL_SetRenderTarget(renderer, nullptr);
    }
    void generate_alpha_mask(SDL_Renderer* renderer)
    {
        SDL_RenderClear(renderer);
        SDL_SetRenderTarget(renderer, alpha_mask);
        SDL_Rect src;
        src.x = 0;
        src.y = 0;
        src.w = dot->width;
        src.h = dot->height;
        SDL_Rect dst;
        dst.w = dot->width;
        dst.h = dot->height;
        for (int dot_y = 0; dot_y < height; ++dot_y)
        {
            for (int dot_x = 0; dot_x < width; ++dot_x)
            {
                dst.x = dot->width * dot_x;
                dst.y = dot->height * dot_y;
                SDL_RenderCopy(renderer, dot->texture, &src, &dst);
            }
        }
    }
    void show_number(SDL_Renderer* renderer, Font* font, int number, int x, int y, TextAlignment align = ALIGN_TOP_LEFT, int max_width = 0, int max_height = 0, bool comma_separated = false, bool trailing_zero = false, int decimal_place = 0)
    {
        std::string number_string = std::to_string(number);
        if (trailing_zero)
        {
            number_string.append("0");
        }
        if (decimal_place)
        {
            number_string.insert(number_string.end() - decimal_place, '.');
        }
        if (comma_separated)
        {
            for (int c = number_string.size() - decimal_place - !!decimal_place - 3; c > 0; c -= 3)
            {
                number_string.insert(number_string.begin() + c, ',');
            }
        }
        show_text(renderer, font, number_string, x, y, align, max_width, max_height);
    }
    void show_number(SDL_Renderer* renderer, Font* font, double number, int x, int y, TextAlignment align = ALIGN_TOP_LEFT, int max_width = 0, int max_height = 0, bool comma_separated = false, int precision = -1)
    {
        std::string number_string = std::to_string(number);
        int decimal_character = 0;
        for (unsigned int c = 0; c < number_string.size(); ++c)
        {
            if (number_string[c] == '.')
            {
                decimal_character = c;
                break;
            }
        }
        if (precision != -1)
        {
            number_string.resize(decimal_character + precision + !!precision);
        }
        if (comma_separated)
        {
            for (int c = decimal_character - 3; c > 0; c -= 3)
            {
                number_string.insert(number_string.begin() + c, ',');
            }
        }
        show_text(renderer, font, number_string, x, y, align, max_width, max_height);
    }
    void show_text(SDL_Renderer* renderer, FontSize* font, std::string text, int x, int y, TextAlignment align = ALIGN_TOP_LEFT)
    {
        SDL_SetRenderTarget(renderer, dot_matrix);
        font->show_text(renderer, text, x, y, align);
        SDL_SetRenderTarget(renderer, nullptr);
    }
    void show_text(SDL_Renderer* renderer, Font* font, std::string text, int x, int y, TextAlignment align = ALIGN_TOP_LEFT, int max_width = 0, int max_height = 0)
    {
        SDL_SetRenderTarget(renderer, dot_matrix);
        font->show_text(renderer, text, x, y, align, max_width, max_height);
        SDL_SetRenderTarget(renderer, nullptr);
    }
    void show_sprite(SDL_Renderer* renderer, Sprite* sprite, double x, double y, double angle, int animation, int frame)
    {
        SDL_SetRenderTarget(renderer, dot_matrix);
        sprite->dmd_render(renderer, x, y, angle, animation, frame);
        SDL_SetRenderTarget(renderer, nullptr);
    }
    void show_point(SDL_Renderer* renderer, int x, int y, int brightness = 255, int alpha = 255)
    {
        SDL_SetRenderTarget(renderer, dot_matrix);
        SDL_SetRenderDrawColor(renderer, brightness, brightness, brightness, alpha);
        SDL_RenderDrawPoint(renderer, x, y);
        SDL_SetRenderTarget(renderer, nullptr);
    }
    void show_line(SDL_Renderer* renderer, int x1, int y1, int x2, int y2, int brightness = 255, int alpha = 255)
    {
        SDL_SetRenderTarget(renderer, dot_matrix);
        SDL_SetRenderDrawColor(renderer, brightness, brightness, brightness, alpha);
        SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
        SDL_SetRenderTarget(renderer, nullptr);
    }
    SpriteInstance* add_sprite(Sprite* sprite, int x, int y)
    {
        //TODO (Phase X) test
        SpriteInstance* instance = new SpriteInstance(sprite, x, y, 0.0);
        sprite_instances.push_back(instance);
        return instance;
    }
    void sprites_present(SDL_Renderer* renderer, double ms)
    {
        //TODO (Phase X) test
        SDL_SetRenderTarget(renderer, dot_matrix);
        for (unsigned int i = 0; i < sprite_instances.size(); ++i)
        {
            sprite_instances[i]->dmd_render(renderer);
            sprite_instances[i]->pass_time(ms);
        }
        SDL_SetRenderTarget(renderer, nullptr);
    }
    void render(SDL_Renderer* renderer, double x, double y, double scale = 1.0, int alpha = 255)
    {
        if (x + width * dot->width * scale < 0 || x > DISPLAY_WIDTH || y + height * dot->height * scale < 0 || y > DISPLAY_HEIGHT)
        {
            return;
        }
        SDL_Rect src;
        src.x = 0;
        src.y = 0;
        src.w = width;
        src.h = height;
        SDL_Rect dst;
        dst.x = 0;
        dst.y = 0;
        dst.w = width * dot->width;
        dst.h = height * dot->height;
        SDL_SetRenderTarget(renderer, intermediate);
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);
        if (invert)
        {
            SDL_SetTextureBlendMode(dot_matrix, SDL_ComposeCustomBlendMode(SDL_BLENDFACTOR_ONE, SDL_BLENDFACTOR_ONE, SDL_BLENDOPERATION_REV_SUBTRACT, SDL_BLENDFACTOR_ONE, SDL_BLENDFACTOR_ZERO, SDL_BLENDOPERATION_ADD));
        } else
        {
            SDL_SetTextureBlendMode(dot_matrix, SDL_BLENDMODE_BLEND);
        }
        SDL_RenderCopy(renderer, dot_matrix, &src, &dst);
        SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, unlit);
        SDL_RenderFillRect(renderer, nullptr);
        SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_NONE);

        SDL_SetTextureColorMod(intermediate, lit_r, lit_g, lit_b);
        SDL_SetTextureAlphaMod(intermediate, alpha);
        src = dst;
        SDL_FRect dstf;
        dstf.x = x;
        dstf.y = y;
        dstf.w = dst.w * scale;
        dstf.h = dst.h * scale;

        SDL_RenderCopy(renderer, alpha_mask, nullptr, nullptr);

        SDL_SetRenderTarget(renderer, nullptr);
        SDL_RenderCopyExF(renderer, intermediate, &src, &dstf, 0.0, nullptr, SDL_FLIP_NONE);
    }
};

enum AlphanumericSpacing
{
    SPACING_PREFIX,
    SPACING_STANDARD,
    SPACING_SUFFIX
};

class AlphanumericDisplay
{
public:
    Sprite* segments = nullptr;
    int lit_r = 255;
    int lit_g = 128;
    int lit_b = 0;
    int unlit = 16;
    int segment_count;
    int character_count;
    int character_width;
    int* segments_lit;
    uint32_t protocol[128] = {0}; //bits correspond to lit segments starting from least significant bit
    SDL_Texture* intermediate = nullptr;
    AlphanumericSpacing spacing_protocol[128] = {SPACING_STANDARD};
    AlphanumericDisplay(SDL_Renderer* renderer, const char filename[], const char protocol_filename[], int width, int height, int segment_count, int character_count) : segment_count(segment_count), character_count(character_count), character_width(width)
    {
        segments = new Sprite(renderer, filename, 0, 0, width, height);
        segments_lit = new int[character_count];
        intermediate = SDL_CreateTexture(renderer, 0, SDL_TEXTUREACCESS_TARGET, character_count * width, height);
        SDL_SetTextureBlendMode(intermediate, SDL_BLENDMODE_BLEND);
        std::ifstream protocol_data(protocol_filename);
        std::string line;
        std::string word;
        std::getline(protocol_data, line);
        std::unordered_map<std::string, int> segment_names;
        int s = 0;
        std::istringstream linestream(line);
        while (linestream >> word)
        {
            segment_names.insert({word, s});
            ++s;
        }
        int i = 0;
        while (i < 128)
        {
            std::getline(protocol_data, line);
            std::istringstream linestream(line);
            spacing_protocol[i] = SPACING_STANDARD;
            while (linestream >> word)
            {
                if (segment_names.find(word) != segment_names.end())
                {
                    protocol[i] |= 1 << segment_names[word];
                } else if (word == "PREFIX")
                {
                    spacing_protocol[i] = SPACING_PREFIX;
                } else if (word == "SUFFIX")
                {
                    spacing_protocol[i] = SPACING_SUFFIX;
                }
            }
            ++i;
        }
        clear_display();
    }
    ~AlphanumericDisplay()
    {
        delete segments;
        delete segments_lit;
        SDL_DestroyTexture(intermediate);
    }
    int get_width(std::string text)
    {
        int width = 0;
        bool on_prefix = false;
        for (unsigned int i = 0; i < text.size(); ++i)
        {
            int c = text[i];
            switch (spacing_protocol[c])
            {
            case SPACING_PREFIX:
                if (!on_prefix)
                {
                    ++width;
                }
                on_prefix = true;
                break;
            case SPACING_STANDARD:
                if (!on_prefix)
                {
                    ++width;
                }
                on_prefix = false;
                break;
            case SPACING_SUFFIX:
                if (width == 0)
                {
                    ++width;
                }
                on_prefix = false;
                break;
            }
        }
        return width;
    }
    void clear_display(bool setting = false)
    {
        for (int i = 0; i < character_count; ++i)
        {
            segments_lit[i] = 0;
        }
    }
    void show_character(int position, int character, bool erase_existing = false)
    {
        if (position >= 0 && position < character_count)
        {
            if (erase_existing)
            {
                segments_lit[position] = 0;
            }
            segments_lit[position] |= protocol[character];
        }
    }
    void show_number(int number, int position, TextAlignment align = ALIGN_LEFT, int max_width = 0, bool comma_separated = false, bool trailing_zero = false, int decimal_place = 0)
    {
        std::string number_string = std::to_string(number);
        if (trailing_zero)
        {
            number_string.append("0");
        }
        if (decimal_place)
        {
            number_string.insert(number_string.end() - decimal_place, '.');
        }
        if (comma_separated)
        {
            for (int c = number_string.size() - decimal_place - !!decimal_place - 3; c > 0; c -= 3)
            {
                number_string.insert(number_string.begin() + c, ',');
            }
        }
        show_text(number_string, position, align, max_width);
    }
    void show_number(double number, int position, TextAlignment align = ALIGN_LEFT, int max_width = 0, bool comma_separated = false, int precision = -1)
    {
        std::string number_string = std::to_string(number);
        int decimal_character = 0;
        for (unsigned int c = 0; c < number_string.size(); ++c)
        {
            if (number_string[c] == '.')
            {
                decimal_character = c;
                break;
            }
        }
        if (precision != -1)
        {
            number_string.resize(decimal_character + precision + !!precision);
        }
        if (comma_separated)
        {
            for (int c = decimal_character - 3; c > 0; c -= 3)
            {
                number_string.insert(number_string.begin() + c, ',');
            }
        }
        show_text(number_string, position, align, max_width);
    }
    void show_text(std::string text, int position, TextAlignment align = ALIGN_LEFT, int max_width = 0)
    {
        if (align == ALIGN_RIGHT || align == ALIGN_TOP_RIGHT || align == ALIGN_BOTTOM_RIGHT)
        {
            position -= get_width(text);
        } else if (align == ALIGN_CENTER || align == ALIGN_TOP || align == ALIGN_BOTTOM)
        {
            position -= get_width(text) / 2;
        }
        for (unsigned int i = 0; i < text.size(); ++i)
        {
            if (max_width && static_cast<int>(i) >= max_width)
            {
                break;
            }
            int c = text[i];
            AlphanumericSpacing c_spacing = spacing_protocol[c];
            if (c_spacing == SPACING_SUFFIX)
            {
                --position;
            }
            show_character(position, c);
            if (c_spacing != SPACING_PREFIX)
            {
                ++position;
            }
        }
    }
    void set_segment(int position, int segment, bool setting = true)
    {
        if (position < 0 && position <= character_count)
        {
            int mask = 1 << segment;
            if (setting)
            {
                segments_lit[position] |= mask;
            } else
            {
                segments_lit[position] &= ~mask;
            }
        }
    }
    void render(SDL_Renderer* renderer, double x, double y, double scale = 1.0, int alpha = 255)
    {
        SDL_SetRenderTarget(renderer, intermediate);
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
        SDL_RenderClear(renderer);
        for (int i = 0; i < character_count; ++i)
        {
            for (int j = 0; j < segment_count; ++j)
            {
                if (1 << j & segments_lit[i])
                {
                    SDL_SetTextureColorMod(segments->texture, lit_r, lit_g, lit_b);
                } else
                {
                    SDL_SetTextureColorMod(segments->texture, unlit, unlit, unlit);
                }
                segments->dmd_render(renderer, character_width * i, 0.0, 0.0, 0, j);
            }
        }

        SDL_SetRenderTarget(renderer, nullptr);
        SDL_Rect src;
        src.x = 0;
        src.y = 0;
        src.w = character_count * character_width;
        src.h = segments->height;
        SDL_FRect dst;
        dst.x = x;
        dst.y = y;
        dst.w = src.w * scale;
        dst.h = src.h * scale;
        SDL_SetTextureColorMod(intermediate, lit_r, lit_g, lit_b);
        SDL_SetTextureAlphaMod(intermediate, alpha);
        SDL_RenderCopyF(renderer, intermediate, &src, &dst);
    }
};

enum Event
{
    LOAD,
    EVERY_TICK,
    EVERY_FRAME,
    END,
    FUNCTION,
    INPUT_DOWN,
    INPUT_UP,
    COLLISION,
    SWITCH_DOWN,
    SWITCH_UP,
    TILT_BOB_CONTACT,
    SPIN_CW,
    SPIN_CCW,
    TIMER_DONE,
    REEL_RESET_DONE,
    REEL_ROLLOVER,
};

class ScoreReel
{
public:
    int id = -1;
    int reel_count;
    int width;
    int height;
    int base;
    double roll_time; //seconds
    SDL_Texture* texture;
    int* showing;
    int* to_show;
    double* time_showing; //seconds
    bool resetting = false;
    Sound* roll_sound = nullptr;
    bool invert_roll = false;
    bool lit = false;
    SDL_Texture* shading_texture = nullptr;
    std::deque<std::pair<Event, int>> trigger_queue;
    ScoreReel(SDL_Renderer* renderer, int reel_count, const char filename[], int width, int height, int base = 10, double roll_time = 0.120) : reel_count(reel_count), width(width), height(height), base(base), roll_time(roll_time)
    {
        SDL_Surface* surface = IMG_Load(filename);
        texture = SDL_CreateTextureFromSurface(renderer, surface);
        SDL_FreeSurface(surface);
        showing = new int[reel_count];
        to_show = new int[reel_count];
        time_showing = new double[reel_count];
        for (int i = 0; i < reel_count; ++i)
        {
            showing[i] = 0;
            to_show[i] = 0;
            time_showing[i] = roll_time;
        }
    }
    ~ScoreReel()
    {
        delete[] showing;
        delete[] time_showing;
        SDL_DestroyTexture(texture);
    }
    void show_number(int number)
    {
        if (number == 0)
        {
            resetting = true;
        }
        int reel = reel_count;
        while(reel--)
        {
            to_show[reel] = number % base;
            number /= base;
        }
    }
    void pass_time(double ms)
    {
        bool stopped = true;
        for (int i = 0; i < reel_count; ++i)
        {
            if (showing[i] != to_show[i] && time_showing[i] >= roll_time && (resetting || i == reel_count - 1 || to_show[i + 1] >= showing[i + 1]))
            {
                if (i == 0 && (showing[i] + 1) == base && !resetting)
                {
                    trigger_queue.push_back(std::pair<Event, int>(REEL_ROLLOVER, id));
                }
                showing[i] = (showing[i] + 1) % base;
                time_showing[i] = 0.0;
                if (roll_sound)
                {
                    roll_sound->play();
                }
            }
            if (time_showing[i] < roll_time)
            {
                stopped = false;
            }
            time_showing[i] += ms / 1000.0;
        }
        if (stopped && resetting)
        {
            trigger_queue.push_back(std::pair<Event, int>(REEL_RESET_DONE, id));
            resetting = false;
        }
    }
    void add_sound(const char filename[])
    {
        roll_sound = new Sound(filename, TRACK_PHYSICAL);
    }
    void add_shading_texture(SDL_Renderer* renderer, const char filename[])
    {
        if (shading_texture)
        {
            SDL_DestroyTexture(shading_texture);
        }
        SDL_Surface* surface = IMG_Load(filename);
        shading_texture = SDL_CreateTextureFromSurface(renderer, surface);
        SDL_FreeSurface(surface);
        SDL_SetTextureBlendMode(shading_texture, SDL_BLENDMODE_BLEND);
    }
    void render(SDL_Renderer* renderer, double x, double y, double scale = 1.0, int alpha = 255)
    {
        for (int i = 0; i < reel_count; ++i)
        {
            SDL_FRect dst;
            dst.x = x + width * scale * i;
            dst.y = y;
            dst.w = width * scale;
            dst.h = height * scale;
            SDL_Rect src;
            src.x = 0;
            if (time_showing[i] >= roll_time)
            {
                src.y = height * (base - showing[i]);
            } else
            {
                src.y = height * (base - (showing[i] + base - 1) % base) - (time_showing[i] * time_showing[i] / roll_time / roll_time) * height;
            }
            src.w = width;
            src.h = height;
            if (invert_roll)
            {
                src.y = height * base - src.y;
            }
            SDL_SetTextureAlphaMod(texture, alpha);
            SDL_RenderCopyF(renderer, texture, &src, &dst);
            if (shading_texture)
            {
                src.y = 0;
                if (lit)
                {
                    src.x = src.w;
                } else
                {
                    src.x = 0;
                }
                SDL_RenderCopyF(renderer, shading_texture, &src, &dst);
            }
        }
    }
};

class Backglass
{
public:
    Point top_left = Point(-13.5 + 10.125, -27.0);
    Point dimensions = Point(27.0, 27.0);
    std::vector<SpriteInstance*> sprites;
    std::vector<std::pair<DotMatrixDisplay*, Point>> dot_matrix_displays;
    std::vector<std::pair<AlphanumericDisplay*, Point>> alphanumeric_displays;
    std::vector<std::pair<ScoreReel*, Point>> score_reels;
    std::deque<std::pair<Event, int>> trigger_queue;
    Backglass()
    {

    }
    ~Backglass()
    {
        for (unsigned int i = 0; i < sprites.size(); ++i)
        {
            delete sprites[i];
        }
    }
    void add_element(Sprite* sprite, double x, double y, double angle)
    {
        sprites.push_back(new SpriteInstance(sprite, top_left.x + x, top_left.y + y, angle));
    }
    void add_element(DotMatrixDisplay* dot_matrix_display, double x, double y)
    {
        dot_matrix_displays.push_back(std::pair<DotMatrixDisplay*, Point>(dot_matrix_display, Point(top_left.x + x, top_left.y + y)));
    }
    void add_element(AlphanumericDisplay* alphanumeric_display, double x, double y)
    {
        alphanumeric_displays.push_back(std::pair<AlphanumericDisplay*, Point>(alphanumeric_display, Point(top_left.x + x, top_left.y + y)));
    }
    void add_element(ScoreReel* score_reel, double x, double y)
    {
        score_reels.push_back(std::pair<ScoreReel*, Point>(score_reel, Point(top_left.x + x, top_left.y + y)));
    }
    void pass_time(double ms)
    {
        for (unsigned int i = 0; i < sprites.size(); ++i)
        {
            sprites[i]->pass_time(ms);
        }
        for (unsigned int i = 0; i < score_reels.size(); ++i)
        {
            score_reels[i].first->pass_time(ms);
            pump_triggers(score_reels[i].first);
        }
    }
    void pump_triggers(ScoreReel* reel)
    {
        while (reel->trigger_queue.size())
        {
            trigger_queue.push_back(reel->trigger_queue.front());
            reel->trigger_queue.pop_front();
        }
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        if (view->debug)
        {
            SDL_Rect boundary;
            Point boundary_top_left = Point(0.0, 0.0);
            view->to_screen_space(&top_left, &boundary_top_left);
            boundary.x = boundary_top_left.x;
            boundary.y = boundary_top_left.y;
            boundary.w = dimensions.x * view->zoom * PIXELS_PER_INCH;
            boundary.h = dimensions.y * view->zoom * PIXELS_PER_INCH;
            SDL_SetRenderDrawColor(renderer, 128, 128, 128, 255);
            SDL_RenderDrawRect(renderer, &boundary);
        }
        for (unsigned int i = 0; i < sprites.size(); ++i)
        {
            sprites[i]->render(renderer, view);
        }
        for (unsigned int i = 0; i < dot_matrix_displays.size(); ++i)
        {
            Point screen_position = Point(0.0, 0.0);
            dot_matrix_displays[i].first->sprites_present(renderer, MS_PER_FRAME);
            view->to_screen_space(&dot_matrix_displays[i].second, &screen_position);
            dot_matrix_displays[i].first->render(renderer, screen_position.x, screen_position.y, view->zoom, 255);
        }
        for (unsigned int i = 0; i < alphanumeric_displays.size(); ++i)
        {
            Point screen_position = Point(0.0, 0.0);
            view->to_screen_space(&alphanumeric_displays[i].second, &screen_position);
            alphanumeric_displays[i].first->render(renderer, screen_position.x, screen_position.y, view->zoom, 255);
        }
        for (unsigned int i = 0; i < score_reels.size(); ++i)
        {
            Point screen_position = Point(0.0, 0.0);
            view->to_screen_space(&score_reels[i].second, &screen_position);
            score_reels[i].first->render(renderer, screen_position.x, screen_position.y, view->zoom, 255);
        }
    }
    void HUD_render(SDL_Renderer* renderer, View* view)
    {
        //TODO margins
        if (!view->HUD)
        {
            return;
        }
        double x = 0.0;
        double y = 0.0;
        double scale = HUD_SCALE;
        double alpha = 255;
        for (unsigned int i = 0; i < dot_matrix_displays.size(); ++i)
        {
            DotMatrixDisplay* DMD = dot_matrix_displays[i].first;
            DMD->render(renderer, x, y, scale, alpha);
            y += DMD->height * DMD->dot->height * scale;
        }
        for (unsigned int i = 0; i < alphanumeric_displays.size(); ++i)
        {
            AlphanumericDisplay* AND = alphanumeric_displays[i].first;
            AND->render(renderer, x, y, scale, alpha);
            y += AND->segments->height * scale;
        }
        for (unsigned int i = 0; i < score_reels.size(); ++i)
        {
            ScoreReel* reel = score_reels[i].first;
            reel->render(renderer, x, y, scale, alpha);
            y += reel->height * scale;
        }
    }
};

enum Command
{
    //View
    DEBUG_RENDERER,
    SPRITE_RENDERER,
    HUD_RENDERER,
    SET_VIEW,
    GET_VIEW,
    //Design
    //TODO turtle commands
    TABLE,
    LAYER,
    POINT,
    LINE,
    ARC,
    CIRCLE,
    TURTLE_TELEPORT,
    TURTLE_SET_HEADING,
    TURTLE_FORWARD,
    TURTLE_TURN,
    TURTLE_POINT,
    TURTLE_LINE,
    TURTLE_ARC,
    TURTLE_GET_POSITION,
    BALL,
    STANDARD_BALL,
    MAGNET,
    SPRITE,
    DMD_SPRITE, //Note: no option to remove DMD sprites is implemented, as it is deemed sufficient to set them as invisible
    BACKGROUND,
    SPINNER,
    TILT_BOB,
    //Manipulate
    GRADIENT_VECTOR,
    BALL_RETURN,
    COMPOSITE,
    MOBILE,
    SOLID,
    WOOD,
    METAL,
    PLASTIC,
    RUBBER,
    RUBBER_BAND,
    WIRE,
    SWITCH,
    GATE,
    MAGNET_POWER,
    SYNC_MAGNET,
    LAYER_PORTAL,
    TANGIBLE,
    LIVE_CATCH,
    LIVE_CATCH_PROPERTIES,
    SET_SENSITIVITY,
    SET_SPRITE_VISIBILITY,
    TRANSLATE_SPRITE,
    ROTATE_SPRITE,
    SET_BALL_SPRITE_ROTATION,
    SYNC_SPRITE,
    SET_SPRITE_ANIMATION,
    SET_SPRITE_FRAME,
    LOOP_SPRITE_ANIMATION,
    LOOP_SYNCHRONOUS_SPRITE_ANIMATION,
    RESET_TILT_BOB,
    SET_DETECTION_NORMAL,
    //Action
    NUDGE,
    HALT_NUDGE,
    BUMPER_ACTION,
    SLINGSHOT_ACTION,
    KICKOUT_ACTION,
    SET_MOBILE_VELOCITY,
    SET_MOBILE_ACCELERATION,
    SET_MOBILE_ANGULAR_VELOCITY,
    SET_MOBILE_ANGULAR_ACCELERATION,
    SET_MOBILE_SPRING_CONSTANT,
    SET_MOBILE_ANGULAR_SPRING_CONSTANT,
    MOBILE_BOUNCE,
    GET_MOBILE_POSITION,
    GET_MOBILE_VELOCITY,
    GET_MOBILE_ACCELERATION,
    GET_MOBILE_ANGLE,
    GET_MOBILE_ANGULAR_VELOCITY,
    GET_MOBILE_ANGULAR_ACCELERATION,
    //God
    SET_BALL_POSITION,
    SET_BALL_VELOCITY,
    SET_BALL_LAYER,
    SET_BALL_TANGIBILITY,
    GET_BALL_POSITION,
    GET_BALL_VELOCITY,
    GET_BALL_LAYER,
    GET_BALL_TANGIBILITY,
    //Input
    GET_DIGITAL_INPUT,
    GET_ANALOG_INPUT,
    //HUD
    DOT_MATRIX_DISPLAY,
    ALPHANUMERIC_DISPLAY,
    SCORE_REEL,
    DMD_SHOW_TEXT,
    DMD_SHOW_NUMBER,
    DMD_SHOW_SPRITE,
    DMD_SHOW_LINE,
    DMD_SHOW_POINT,
    DMD_INVERT,
    DMD_CLEAR,
    AND_SHOW_TEXT,
    AND_SHOW_NUMBER,
    AND_SET_SEGMENT,
    AND_CLEAR,
    REEL_SHOW_NUMBER,
    REEL_RESET,
    REEL_ADD_SHADING,
    REEL_LIGHT,
    //Font
    LOAD_FONT_SIZE,
    LOAD_MONOSPACED_FONT_SIZE,
    CREATE_FONT,
    ADD_FONT_SIZE,
    //Backglass
    BACKGLASS,
    ADD_DOT_MATRIX_DISPLAY,
    ADD_ALPHANUMERIC_DISPLAY,
    ADD_SCORE_REEL,
    DESIGNATE_BACKGLASS,
    EMBED_BACKGLASS,
    //Sound
    PLAY_SOUND,
    STOP_SOUND,
    PLAY_MUSIC,
    TRANSITION_MUSIC,
    STOP_MUSIC,
    //Variable
    INTEGER,
    DOUBLE,
    STRING,
    ARRAY,
    SET_INTEGER,
    SET_DOUBLE,
    SET_STRING,
    SET_ARRAY_ELEMENT,
    GET_ARRAY_ELEMENT,
    STRING_LENGTH,
    STRING_CHARACTER,
    STRING_APPEND,
    STRING_CONCATENATE,
    RANDOM,
    ADD,
    SUBTRACT,
    MULTIPLY,
    DIVIDE,
    MODULO,
    D_ADD,
    D_SUBTRACT,
    D_MULTIPLY,
    D_DIVIDE,
    AND,
    OR,
    XOR,
    NOT,
    BITWISE_AND,
    BITWISE_OR,
    BITWISE_XOR,
    BITWISE_NOT,
    GREATER,
    LESS,
    EQUAL,
    D_GREATER,
    D_LESS,
    D_EQUAL,
    SHIFT_ARRAY,
    CLEAR_ARRAY,
    INTEGER_TO_DOUBLE,
    DOUBLE_TO_INTEGER,
    INTEGER_TO_STRING,
    DOUBLE_TO_STRING,
    //Timer
    TIMER,
    START_TIMER,
    PAUSE_TIMER,
    RESUME_TIMER,
    TIMER_TO_DOUBLE,
    //Flow control
    IF,
    ELSE,
    ENDIF,
    WHILE,
    WEND,
    BREAK,
    GOSUB,
    RETURN,
    LABEL,
    //File
    LOAD_SPRITESHEET,
    LOAD_SOUND,
    LOAD_MUSIC,
    LOAD_RECORD,
    LOAD_CONFIG,
    LOAD_CONTROLS,
    SAVE_RECORD,
    SAVE_CONFIG,
    SAVE_CONTROLS,
    GET_RECORD,
    QUALIFIES,
    DEFAULT_CONFIG,
    DEFAULT_CONTROLS,
    //Settings
    BIND_KEY,
    BIND_BUTTON,
    BIND_AXIS,
    UNBIND_KEY,
    UNBIND_BUTTON,
    UNBIND_AXIS,
    UNBIND_ALL,
    DISPLAY_MODE,
    RESOLUTION,
    FRAMERATE,
    FLASHING_LIGHTS,
    FLASHING_LIGHTS_ENABLED,
    MUTE,
    MASTER_VOLUME,
    MUSIC_VOLUME,
    SOUND_VOLUME,
    SCALE_UI,
    SCALE_HUD,
    //Console
    PRINT,
    //Macro
    FLIPPER,
    //TODO (Phase VII) SLINGSHOT (should take 3 points and a bevel radius and output ID of active face
};

enum DataType
{
    TYPE_INTEGER,
    TYPE_DOUBLE,
    TYPE_STRING,
    TYPE_INTEGER_VARIABLE,
    TYPE_DOUBLE_VARIABLE,
    TYPE_STRING_VARIABLE,
    TYPE_ARRAY,
    TYPE_NULL,
};

enum Input
{
    INPUT_START,
    INPUT_COIN,
    INPUT_LEFT_FLIPPER,
    INPUT_RIGHT_FLIPPER,
    INPUT_LEFT_UPPER_FLIPPER,
    INPUT_RIGHT_UPPER_FLIPPER,
    INPUT_LEFT_SECONDARY,
    INPUT_RIGHT_SECONDARY,
    INPUT_BALL_ELEVATOR,
    INPUT_PLUNGER,
    INPUT_NUDGE,
    INPUT_MENU_UP,
    INPUT_MENU_DOWN,
    INPUT_MENU_LEFT,
    INPUT_MENU_RIGHT,
    INPUT_MENU_CONFIRM,
    INPUT_MENU_CANCEL,
    INPUT_CAMERA_ZOOM_IN,
    INPUT_CAMERA_ZOOM_OUT,
    INPUT_CAMERA_HOME,
    ANALOG_PLUNGER,
    ANALOG_NUDGE_X,
    ANALOG_NUDGE_Y,
    ANALOG_CAMERA_X,
    ANALOG_CAMERA_Y,
    ANALOG_CAMERA_ZOOM,
    INPUT_COUNT
};

enum GradientMode
{
    GRADIENT_CONSTANT,
    GRADIENT_POLAR,
};

class Syntax
{
public:
    std::unordered_map<std::string, Command> command_table;
    std::unordered_map<Command, std::vector<DataType>*> arg_table;
    std::unordered_map<std::string, Event> event_table;
    std::unordered_map<std::string, int> enum_table;
    Syntax()
    {
        //[SYNTAX LIST]
        command_table["DEBUG_RENDERER"] = DEBUG_RENDERER;
        define_syntax(DEBUG_RENDERER, 1, TYPE_INTEGER);
        command_table["SPRITE_RENDERER"] = SPRITE_RENDERER;
        define_syntax(SPRITE_RENDERER, 1, TYPE_INTEGER);
        command_table["HUD_RENDERER"] = HUD_RENDERER;
        define_syntax(HUD_RENDERER, 1, TYPE_INTEGER);
        command_table["SET_VIEW"] = SET_VIEW;
        define_syntax(SET_VIEW, 4, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE); //x, y, zoom, transition time
        command_table["GET_VIEW"] = GET_VIEW;
        define_syntax(GET_VIEW, 3, TYPE_DOUBLE_VARIABLE, TYPE_DOUBLE_VARIABLE, TYPE_DOUBLE_VARIABLE); //x out, y out, zoom out
        command_table["TABLE"] = TABLE;
        define_syntax(TABLE, 3, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER); //width, height, nudge mode
        command_table["LAYER"] = LAYER;
        define_syntax(LAYER, 10, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //top left x, top left y, width, height, gradient mode, gradient vector x, gradient vector y, center x, center y, ID out
        command_table["POINT"] = POINT;
        define_syntax(POINT, 3, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //x, y, ID out
        command_table["LINE"] = LINE;
        define_syntax(LINE, 5, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //x1, y1, x2, y2, ID out
        command_table["ARC"] = ARC;
        define_syntax(ARC, 6, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //x, y, radius, start angle, end angle, ID out
        command_table["CIRCLE"] = CIRCLE;
        define_syntax(CIRCLE, 4, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //x, y, radius, ID out
        command_table["TURTLE_TELEPORT"] = TURTLE_TELEPORT;
        define_syntax(TURTLE_TELEPORT, 2, TYPE_DOUBLE, TYPE_DOUBLE); //x, y
        command_table["TURTLE_SET_HEADING"] = TURTLE_SET_HEADING;
        define_syntax(TURTLE_SET_HEADING, 1, TYPE_DOUBLE); //heading
        command_table["TURTLE_FORWARD"] = TURTLE_FORWARD;
        define_syntax(TURTLE_FORWARD, 1, TYPE_DOUBLE); //distance
        command_table["TURTLE_TURN"] = TURTLE_TURN;
        define_syntax(TURTLE_TURN, 1, TYPE_DOUBLE); //angle
        command_table["TURTLE_POINT"] = TURTLE_POINT;
        define_syntax(TURTLE_POINT, 1, TYPE_INTEGER_VARIABLE); //ID out
        command_table["TURTLE_LINE"] = TURTLE_LINE;
        define_syntax(TURTLE_LINE, 2, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //length, ID out
        command_table["TURTLE_ARC"] = TURTLE_ARC;
        define_syntax(TURTLE_ARC, 3, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //radius, angle, ID out
        command_table["TURTLE_GET_POSITION"] = TURTLE_GET_POSITION;
        define_syntax(TURTLE_GET_POSITION, 2, TYPE_DOUBLE_VARIABLE, TYPE_DOUBLE_VARIABLE); //x out, y out
        command_table["BALL"] = BALL;
        define_syntax(BALL, 12, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //layer, x, y, radius, mass, moment ratio, magnetic, show spin, spritesheet, contact sound, roll sound, ID out
        command_table["STANDARD_BALL"] = STANDARD_BALL;
        define_syntax(STANDARD_BALL, 4, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //layer, x, y, ID out
        command_table["MAGNET"] = MAGNET;
        define_syntax(MAGNET, 8, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //layer, x, y, magnetism, radius, depth, hole, ID out
        command_table["SPRITE"] = SPRITE;
        define_syntax(SPRITE, 6, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //ID in, layer, above/below ball, x, y, ID out
        command_table["DMD_SPRITE"] = DMD_SPRITE;
        define_syntax(DMD_SPRITE, 5, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //ID in, DMD ID, x, y, ID out
        command_table["BACKGROUND"] = BACKGROUND;
        define_syntax(BACKGROUND, 1, TYPE_STRING); //filename
        command_table["SPINNER"] = SPINNER;
        define_syntax(SPINNER, 9, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //collider ID in, layer, contact point radius, center mass radius, damping coefficient, gravity, sprite ID in, frames per rotation, ID out
        command_table["TILT_BOB"] = TILT_BOB; //TODO ID out
        define_syntax(TILT_BOB, 6, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE); //x, y, tolerance, damping, bounce coefficient, gravity
        command_table["RESET_TILT_BOB"] = RESET_TILT_BOB;
        define_syntax(RESET_TILT_BOB, 1, TYPE_INTEGER); //tilt bob ID
        command_table["SET_DETECTION_NORMAL"] = SET_DETECTION_NORMAL;
        define_syntax(SET_DETECTION_NORMAL, 2, TYPE_INTEGER, TYPE_DOUBLE); //collider ID, detection normal angle
        command_table["GRADIENT_VECTOR"] = GRADIENT_VECTOR;
        define_syntax(GRADIENT_VECTOR, 3, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //layer, x, y
        command_table["BALL_RETURN"] = BALL_RETURN;
        define_syntax(BALL_RETURN, 3, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //layer, x, y
        command_table["COMPOSITE"] = COMPOSITE; //TODO (Phase VIII) Consider alternate version for joining non-contiguous IDs
        define_syntax(COMPOSITE, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //ID in range start, ID in range end, ID out
        command_table["MOBILE"] = MOBILE;
        define_syntax(MOBILE, 12, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //ID in, origin x, origin y, left range, right range, up range, down range, clockwise range, counterclockwise range, translation bounce, rotation bounce, ID out
        command_table["SOLID"] = SOLID;
        define_syntax(SOLID, 5, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER); //ID in, layer, bounce coefficient, friction coefficient, impact sound
        command_table["WOOD"] = WOOD;
        define_syntax(WOOD, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, layer
        command_table["METAL"] = METAL;
        define_syntax(METAL, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, layer
        command_table["PLASTIC"] = PLASTIC;
        define_syntax(PLASTIC, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, layer
        command_table["RUBBER"] = RUBBER;
        define_syntax(RUBBER, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, layer
        command_table["RUBBER_BAND"] = RUBBER_BAND;
        define_syntax(RUBBER_BAND, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, layer
        command_table["WIRE"] = METAL; //WIRE is an alias for METAL
        command_table["SWITCH"] = SWITCH;
        define_syntax(SWITCH, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //ID in, layer, ID out
        command_table["GATE"] = GATE;
        define_syntax(GATE, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, gate direction
        command_table["MAGNET_POWER"] = MAGNET_POWER;
        define_syntax(MAGNET_POWER, 2, TYPE_INTEGER, TYPE_DOUBLE); //ID in, magnet power
        command_table["SYNC_MAGNET"] = SYNC_MAGNET;
        define_syntax(SYNC_MAGNET, 4, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //magnet ID in, mobile collider ID in, sync point x, sync point y
        command_table["LAYER_PORTAL"] = LAYER_PORTAL;
        define_syntax(LAYER_PORTAL, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, destination layer
        command_table["TANGIBLE"] = TANGIBLE;
        define_syntax(TANGIBLE, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, tangibility mode
        command_table["LIVE_CATCH"] = LIVE_CATCH;
        define_syntax(LIVE_CATCH, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, live catch on/off
        command_table["LIVE_CATCH_PROPERTIES"] = LIVE_CATCH_PROPERTIES;
        define_syntax(LIVE_CATCH_PROPERTIES, 3, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //ID in, live catch elasticity, live catch damping
        command_table["SET_SENSITIVITY"] = SET_SENSITIVITY;
        define_syntax(SET_SENSITIVITY, 2, TYPE_INTEGER, TYPE_DOUBLE); //ID in, sensitivity
        command_table["SET_SPRITE_VISIBILITY"] = SET_SPRITE_VISIBILITY;
        define_syntax(SET_SPRITE_VISIBILITY, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, visibility
        command_table["TRANSLATE_SPRITE"] = TRANSLATE_SPRITE;
        define_syntax(TRANSLATE_SPRITE, 4, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //ID in, x, y
        command_table["ROTATE_SPRITE"] = ROTATE_SPRITE;
        define_syntax(ROTATE_SPRITE, 2, TYPE_INTEGER, TYPE_DOUBLE); //ID in, angle
        command_table["SET_BALL_SPRITE_ROTATION"] = SET_BALL_SPRITE_ROTATION;
        define_syntax(SET_BALL_SPRITE_ROTATION, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, true/false
        command_table["SYNC_SPRITE"] = SYNC_SPRITE;
        define_syntax(SYNC_SPRITE, 2, TYPE_INTEGER, TYPE_INTEGER); //sprite ID in, mobile collider ID in
        command_table["SET_SPRITE_ANIMATION"] = SET_SPRITE_ANIMATION;
        define_syntax(SET_SPRITE_ANIMATION, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, animation
        command_table["SET_SPRITE_FRAME"] = SET_SPRITE_FRAME;
        define_syntax(SET_SPRITE_FRAME, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, frame
        command_table["LOOP_SPRITE_ANIMATION"] = LOOP_SPRITE_ANIMATION;
        define_syntax(LOOP_SPRITE_ANIMATION, 6, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_INTEGER); //ID in, animation, start frame, frame count, frames per second, loop count
        command_table["LOOP_SYNCHRONOUS_SPRITE_ANIMATION"] = LOOP_SYNCHRONOUS_SPRITE_ANIMATION;
        define_syntax(LOOP_SYNCHRONOUS_SPRITE_ANIMATION, 6, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //ID in, animation, start frame, frame count, frames per second, sync offset
        command_table["NUDGE"] = NUDGE;
        define_syntax(NUDGE, 3, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE); //velocity x, velocity y, duration
        command_table["HALT_NUDGE"] = HALT_NUDGE;
        define_syntax(HALT_NUDGE, 0); //no arguments
        command_table["BUMPER_ACTION"] = BUMPER_ACTION;
        define_syntax(BUMPER_ACTION, 4, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //layer, collider ID in, added radius, speed
        command_table["SLINGSHOT_ACTION"] = SLINGSHOT_ACTION;
        define_syntax(SLINGSHOT_ACTION, 5, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE); //layer, collider ID in, kicker range, kicker speed, kicker bias
        command_table["KICKOUT_ACTION"] = KICKOUT_ACTION;
        define_syntax(KICKOUT_ACTION, 6, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //origin layer, x, y, destination layer, impulse x, impulse y
        command_table["SET_MOBILE_VELOCITY"] = SET_MOBILE_VELOCITY;
        define_syntax(SET_MOBILE_VELOCITY, 3, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //ID in, velocity x, velocity y
        command_table["SET_MOBILE_ACCELERATION"] = SET_MOBILE_ACCELERATION;
        define_syntax(SET_MOBILE_ACCELERATION, 3, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //ID in, acceleration x, acceleration y
        command_table["SET_MOBILE_ANGULAR_VELOCITY"] = SET_MOBILE_ANGULAR_VELOCITY;
        define_syntax(SET_MOBILE_ANGULAR_VELOCITY, 2, TYPE_INTEGER, TYPE_DOUBLE); //ID in, angular velocity
        command_table["SET_MOBILE_ANGULAR_ACCELERATION"] = SET_MOBILE_ANGULAR_ACCELERATION;
        define_syntax(SET_MOBILE_ANGULAR_ACCELERATION, 2, TYPE_INTEGER, TYPE_DOUBLE); //ID in, angular acceleration
        command_table["SET_MOBILE_SPRING_CONSTANT"] = SET_MOBILE_SPRING_CONSTANT;
        define_syntax(SET_MOBILE_SPRING_CONSTANT, 2, TYPE_INTEGER, TYPE_DOUBLE); //ID in, spring constant (acceleration per inch of displacement)
        command_table["SET_MOBILE_ANGULAR_SPRING_CONSTANT"] = SET_MOBILE_ANGULAR_SPRING_CONSTANT;
        define_syntax(SET_MOBILE_ANGULAR_SPRING_CONSTANT, 2, TYPE_INTEGER, TYPE_DOUBLE); //ID in, angular spring constant (acceleration per radian of displacement)
        command_table["MOBILE_BOUNCE"] = MOBILE_BOUNCE;
        define_syntax(MOBILE_BOUNCE, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, bounce enabled bool
        command_table["GET_MOBILE_POSITION"] = GET_MOBILE_POSITION;
        define_syntax(GET_MOBILE_POSITION, 3, TYPE_INTEGER, TYPE_DOUBLE_VARIABLE, TYPE_DOUBLE_VARIABLE); //ID in, x out, y out
        command_table["GET_MOBILE_VELOCITY"] = GET_MOBILE_VELOCITY;
        define_syntax(GET_MOBILE_VELOCITY, 3, TYPE_INTEGER, TYPE_DOUBLE_VARIABLE, TYPE_DOUBLE_VARIABLE); //ID in, velocity x out, velocity y out
        command_table["GET_MOBILE_ACCELERATION"] = GET_MOBILE_ACCELERATION;
        define_syntax(GET_MOBILE_ACCELERATION, 3, TYPE_INTEGER, TYPE_DOUBLE_VARIABLE, TYPE_DOUBLE_VARIABLE); //ID in, acceleration x out, acceleration y out
        command_table["GET_MOBILE_ANGLE"] = GET_MOBILE_ANGLE;
        define_syntax(GET_MOBILE_ANGLE, 2, TYPE_INTEGER, TYPE_DOUBLE_VARIABLE); //ID in, angle out
        command_table["GET_MOBILE_ANGULAR_VELOCITY"] = GET_MOBILE_ANGULAR_VELOCITY;
        define_syntax(GET_MOBILE_ANGULAR_VELOCITY, 2, TYPE_INTEGER, TYPE_DOUBLE_VARIABLE); //ID in, angular velocity out
        command_table["GET_MOBILE_ANGULAR_ACCELERATION"] = GET_MOBILE_ANGULAR_ACCELERATION;
        define_syntax(GET_MOBILE_ANGULAR_ACCELERATION, 2, TYPE_INTEGER, TYPE_DOUBLE_VARIABLE); //ID in, angular acceleration out
        command_table["SET_BALL_POSITION"] = SET_BALL_POSITION;
        define_syntax(SET_BALL_POSITION, 3, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //ID in, x, y
        command_table["SET_BALL_VELOCITY"] = SET_BALL_VELOCITY; //TODO (Phase VIII) Rework to include angular velocity
        define_syntax(SET_BALL_VELOCITY, 3, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //ID in, velocity x, velocity y
        command_table["SET_BALL_LAYER"] = SET_BALL_LAYER;
        define_syntax(SET_BALL_LAYER, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, layer
        command_table["SET_BALL_TANGIBILITY"] = SET_BALL_TANGIBILITY;
        define_syntax(SET_BALL_TANGIBILITY, 2, TYPE_INTEGER, TYPE_INTEGER); //ID in, tangibility mode
        command_table["GET_BALL_POSITION"] = GET_BALL_POSITION;
        define_syntax(GET_BALL_POSITION, 3, TYPE_INTEGER, TYPE_DOUBLE_VARIABLE, TYPE_DOUBLE_VARIABLE); //ID in, x out, y out
        command_table["GET_BALL_VELOCITY"] = GET_BALL_VELOCITY; //TODO (Phase VIII) Rework to include angular velocity
        define_syntax(GET_BALL_VELOCITY, 3, TYPE_INTEGER, TYPE_DOUBLE_VARIABLE, TYPE_DOUBLE_VARIABLE); //ID in, velocity x out, velocity y out
        command_table["GET_BALL_LAYER"] = GET_BALL_LAYER;
        define_syntax(GET_BALL_LAYER, 2, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //ID in, layer out
        command_table["GET_BALL_TANGIBILITY"] = GET_BALL_TANGIBILITY;
        define_syntax(GET_BALL_TANGIBILITY, 2, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //ID in, tangibility mode out
        command_table["GET_DIGITAL_INPUT"] = GET_DIGITAL_INPUT;
        define_syntax(GET_DIGITAL_INPUT, 2, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //input, digital value out
        command_table["GET_ANALOG_INPUT"] = GET_ANALOG_INPUT;
        define_syntax(GET_ANALOG_INPUT, 2, TYPE_INTEGER, TYPE_DOUBLE_VARIABLE); //input, analog value out
        command_table["DOT_MATRIX_DISPLAY"] = DOT_MATRIX_DISPLAY;
        define_syntax(DOT_MATRIX_DISPLAY, 7, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //resolution x, resolution y, r, g, b, dot sprite, ID out
        command_table["ALPHANUMERIC_DISPLAY"] = ALPHANUMERIC_DISPLAY;
        define_syntax(ALPHANUMERIC_DISPLAY, 10, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_STRING, TYPE_STRING, TYPE_INTEGER_VARIABLE); //character count, segment count, r, g, b, character width, character height, segments filename, protocol filename, ID out
        command_table["SCORE_REEL"] = SCORE_REEL;
        define_syntax(SCORE_REEL, 9, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_STRING, TYPE_STRING, TYPE_INTEGER, TYPE_DOUBLE, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //reel count, reel width, reel height, filename, roll sound filename, base, roll time, invert roll, ID out
        command_table["DMD_SHOW_TEXT"] = DMD_SHOW_TEXT;
        define_syntax(DMD_SHOW_TEXT, 8, TYPE_INTEGER, TYPE_STRING, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER); //DMD ID in, text, font, x, y, max width, max height, alignment
        command_table["DMD_SHOW_NUMBER"] = DMD_SHOW_NUMBER;
        define_syntax(DMD_SHOW_NUMBER, 11, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER); //DMD ID in, number, font, x, y, max width, max height, alignment, decimal place, comma separation, trailing zero
        command_table["DMD_SHOW_SPRITE"] = DMD_SHOW_SPRITE;
        define_syntax(DMD_SHOW_SPRITE, 7, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER); //DMD ID in, sprite ID in, animation, frame, x, y, alignment
        command_table["DMD_SHOW_LINE"] = DMD_SHOW_LINE;
        define_syntax(DMD_SHOW_LINE, 7, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER); //DMD ID in, x1, y1, x2, y2, brightness, alpha
        command_table["DMD_SHOW_POINT"] = DMD_SHOW_POINT;
        define_syntax(DMD_SHOW_POINT, 5, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER); //DMD ID in, x, y, brightness, alpha
        command_table["DMD_INVERT"] = DMD_INVERT;
        define_syntax(DMD_INVERT, 2, TYPE_INTEGER, TYPE_INTEGER); //DMD ID in, invert on/off
        command_table["DMD_CLEAR"] = DMD_CLEAR;
        define_syntax(DMD_CLEAR, 2, TYPE_INTEGER, TYPE_INTEGER); //DMD ID in, brightness
        command_table["AND_SHOW_TEXT"] = AND_SHOW_TEXT;
        define_syntax(AND_SHOW_TEXT, 4, TYPE_INTEGER, TYPE_STRING, TYPE_INTEGER, TYPE_INTEGER); //AND ID in, text, anchor character, alignment
        command_table["AND_SHOW_NUMBER"] = AND_SHOW_NUMBER;
        define_syntax(AND_SHOW_NUMBER, 7, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER); //AND ID in, number, anchor character, alignment, decimal place, comma separation, trailing zero
        command_table["AND_SET_SEGMENT"] = AND_SET_SEGMENT;
        define_syntax(AND_SET_SEGMENT, 4, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER); //AND ID in, character, segment, on/off
        command_table["AND_CLEAR"] = AND_CLEAR;
        define_syntax(AND_CLEAR, 2, TYPE_INTEGER, TYPE_INTEGER); //AND ID in, on/off
        command_table["REEL_SHOW_NUMBER"] = REEL_SHOW_NUMBER;
        define_syntax(REEL_SHOW_NUMBER, 2, TYPE_INTEGER, TYPE_INTEGER); //reel ID in, number
        command_table["REEL_RESET"] = REEL_RESET;
        define_syntax(REEL_RESET, 1, TYPE_INTEGER); //reel ID in
        command_table["REEL_ADD_SHADING"] = REEL_ADD_SHADING;
        define_syntax(REEL_ADD_SHADING, 2, TYPE_INTEGER, TYPE_STRING); //reel ID in, shading texture filename
        command_table["REEL_LIGHT"] = REEL_LIGHT;
        define_syntax(REEL_LIGHT, 2, TYPE_INTEGER, TYPE_INTEGER); //reel ID in, light on/off
        command_table["LOAD_FONT_SIZE"] = LOAD_FONT_SIZE;
        define_syntax(LOAD_FONT_SIZE, 3, TYPE_STRING, TYPE_STRING, TYPE_INTEGER_VARIABLE); //filename, kerning data filename, ID out
        command_table["LOAD_MONOSPACED_FONT_SIZE"] = LOAD_MONOSPACED_FONT_SIZE;
        define_syntax(LOAD_MONOSPACED_FONT_SIZE, 4, TYPE_STRING, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //filename, width, height, ID out
        command_table["CREATE_FONT"] = CREATE_FONT;
        define_syntax(CREATE_FONT, 1, TYPE_INTEGER_VARIABLE); //ID out
        command_table["ADD_FONT_SIZE"] = ADD_FONT_SIZE;
        define_syntax(ADD_FONT_SIZE, 2, TYPE_INTEGER, TYPE_INTEGER); //font ID in, font size ID in
        command_table["BACKGLASS"] = BACKGLASS;
        define_syntax(BACKGLASS, 5, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //x, y, width, height, ID out
        command_table["ADD_DOT_MATRIX_DISPLAY"] = ADD_DOT_MATRIX_DISPLAY;
        define_syntax(ADD_DOT_MATRIX_DISPLAY, 4, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //backglass ID in, DMD ID in, x, y
        command_table["ADD_ALPHANUMERIC_DISPLAY"] = ADD_ALPHANUMERIC_DISPLAY;
        define_syntax(ADD_ALPHANUMERIC_DISPLAY, 4, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //backglass ID in, AND ID in, x, y
        command_table["ADD_SCORE_REEL"] = ADD_SCORE_REEL;
        define_syntax(ADD_SCORE_REEL, 4, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //backglass ID in, reel ID in, x, y
        command_table["DESIGNATE_BACKGLASS"] = DESIGNATE_BACKGLASS;
        define_syntax(DESIGNATE_BACKGLASS, 1, TYPE_INTEGER); //backglass ID in
        command_table["EMBED_BACKGLASS"] = EMBED_BACKGLASS;
        define_syntax(EMBED_BACKGLASS, 2, TYPE_INTEGER, TYPE_INTEGER); //backglass ID in, layer
        command_table["PLAY_SOUND"] = PLAY_SOUND;
        define_syntax(PLAY_SOUND, 2, TYPE_INTEGER, TYPE_INTEGER); //sound ID in, volume
        command_table["STOP_SOUND"] = STOP_SOUND;
        define_syntax(STOP_SOUND, 0); //no arguments; stops all sounds
        command_table["PLAY_MUSIC"] = PLAY_MUSIC;
        define_syntax(PLAY_MUSIC, 2, TYPE_INTEGER, TYPE_INTEGER); //music ID in, volume
        command_table["TRANSITION_MUSIC"] = TRANSITION_MUSIC;
        define_syntax(TRANSITION_MUSIC, 2, TYPE_INTEGER, TYPE_INTEGER); //music ID in, volume
        command_table["STOP_MUSIC"] = STOP_MUSIC;
        define_syntax(STOP_MUSIC, 0); //no arguments; stops music
        command_table["INTEGER"] = INTEGER;
        define_syntax(INTEGER, 2, TYPE_STRING, TYPE_INTEGER); //variable name, initialized value
        command_table["DOUBLE"] = DOUBLE;
        define_syntax(DOUBLE, 2, TYPE_STRING, TYPE_DOUBLE); //variable name, initialized value
        command_table["STRING"] = STRING;
        define_syntax(STRING, 2, TYPE_STRING, TYPE_STRING); //variable name, initialized value
        command_table["ARRAY"] = ARRAY;
        define_syntax(ARRAY, 3, TYPE_STRING, TYPE_INTEGER, TYPE_INTEGER); //array name, length, initialized value
        command_table["SET_INTEGER"] = SET_INTEGER;
        define_syntax(SET_INTEGER, 2, TYPE_INTEGER_VARIABLE, TYPE_INTEGER); //variable, value
        command_table["SET_DOUBLE"] = SET_DOUBLE;
        define_syntax(SET_DOUBLE, 2, TYPE_DOUBLE_VARIABLE, TYPE_DOUBLE); //variable, value
        command_table["SET_STRING"] = SET_STRING;
        define_syntax(SET_STRING, 2, TYPE_STRING_VARIABLE, TYPE_STRING); //variable, value
        command_table["SET_ARRAY_ELEMENT"] = SET_ARRAY_ELEMENT;
        define_syntax(SET_ARRAY_ELEMENT, 3, TYPE_ARRAY, TYPE_INTEGER, TYPE_INTEGER); //array, index, value
        command_table["GET_ARRAY_ELEMENT"] = GET_ARRAY_ELEMENT;
        define_syntax(GET_ARRAY_ELEMENT, 3, TYPE_ARRAY, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //array, index, value out
        command_table["STRING_LENGTH"] = STRING_LENGTH;
        define_syntax(STRING_LENGTH, 2, TYPE_STRING, TYPE_INTEGER_VARIABLE); //string, length out
        command_table["STRING_CHARACTER"] = STRING_CHARACTER;
        define_syntax(STRING_CHARACTER, 3, TYPE_STRING, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //string, index, character out
        command_table["STRING_APPEND"] = STRING_APPEND;
        define_syntax(STRING_APPEND, 2, TYPE_STRING_VARIABLE, TYPE_INTEGER); //string in, character
        command_table["STRING_CONCATENATE"] = STRING_CONCATENATE;
        define_syntax(STRING_CONCATENATE, 2, TYPE_STRING_VARIABLE, TYPE_STRING); //string in, suffix
        command_table["RANDOM"] = RANDOM;
        define_syntax(RANDOM, 2, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //range (excluded max), random number out
        command_table["ADD"] = ADD;
        define_syntax(ADD, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //addend, addend, sum out
        command_table["SUBTRACT"] = SUBTRACT;
        define_syntax(SUBTRACT, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //minuend, subtrahend, difference out
        command_table["MULTIPLY"] = MULTIPLY;
        define_syntax(MULTIPLY, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //factor, factor, product out
        command_table["DIVIDE"] = DIVIDE;
        define_syntax(DIVIDE, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //dividend, divisor, quotient out
        command_table["MODULO"] = MODULO;
        define_syntax(MODULO, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //modulend, modulator, modulus out
        command_table["D_ADD"] = D_ADD;
        define_syntax(D_ADD, 3, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE_VARIABLE); //addend, addend, sum out
        command_table["D_SUBTRACT"] = D_SUBTRACT;
        define_syntax(D_SUBTRACT, 3, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE_VARIABLE); //minuend, subtrahend, difference out
        command_table["D_MULTIPLY"] = D_MULTIPLY;
        define_syntax(D_MULTIPLY, 3, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE_VARIABLE); //factor, factor, product out
        command_table["D_DIVIDE"] = D_DIVIDE;
        define_syntax(D_DIVIDE, 3, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE_VARIABLE); //dividend, divisor, quotient out
        command_table["AND"] = AND;
        define_syntax(AND, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //A, B, out
        command_table["OR"] = OR;
        define_syntax(OR, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //A, B, out
        command_table["XOR"] = XOR;
        define_syntax(XOR, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //A, B, out
        command_table["NOT"] = NOT;
        define_syntax(NOT, 2, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //in, out
        command_table["BITWISE_AND"] = BITWISE_AND;
        define_syntax(BITWISE_AND, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //A, B, out
        command_table["BITWISE_OR"] = BITWISE_OR;
        define_syntax(BITWISE_OR, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //A, B, out
        command_table["BITWISE_XOR"] = BITWISE_XOR;
        define_syntax(BITWISE_XOR, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //A, B, out
        command_table["BITWISE_NOT"] = BITWISE_NOT;
        define_syntax(BITWISE_NOT, 2, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //in, out
        command_table["GREATER"] = GREATER;
        define_syntax(GREATER, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //A, B, out
        command_table["LESS"] = LESS;
        define_syntax(LESS, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //A, B, out
        command_table["EQUAL"] = EQUAL;
        define_syntax(EQUAL, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //A, B, out
        command_table["D_GREATER"] = D_GREATER;
        define_syntax(D_GREATER, 3, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //A, B, out
        command_table["D_LESS"] = D_LESS;
        define_syntax(D_LESS, 3, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //A, B, out
        command_table["D_EQUAL"] = D_EQUAL;
        define_syntax(D_EQUAL, 3, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //A, B, out
        command_table["SHIFT_ARRAY"] = SHIFT_ARRAY;
        define_syntax(SHIFT_ARRAY, 2, TYPE_ARRAY, TYPE_INTEGER); //array, shift length
        command_table["CLEAR_ARRAY"] = CLEAR_ARRAY;
        define_syntax(CLEAR_ARRAY, 2, TYPE_ARRAY, TYPE_INTEGER); //array, value
        command_table["INTEGER_TO_DOUBLE"] = INTEGER_TO_DOUBLE;
        define_syntax(INTEGER_TO_DOUBLE, 2, TYPE_INTEGER, TYPE_DOUBLE_VARIABLE); //integer, double out
        command_table["DOUBLE_TO_INTEGER"] = DOUBLE_TO_INTEGER;
        define_syntax(DOUBLE_TO_INTEGER, 2, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //double, integer out
        command_table["INTEGER_TO_STRING"] = INTEGER_TO_STRING;
        define_syntax(INTEGER_TO_STRING, 2, TYPE_INTEGER, TYPE_STRING_VARIABLE); //integer, string out
        command_table["DOUBLE_TO_STRING"] = DOUBLE_TO_STRING;
        define_syntax(DOUBLE_TO_STRING, 2, TYPE_DOUBLE, TYPE_STRING_VARIABLE); //double, string out
        command_table["TIMER"] = TIMER;
        define_syntax(TIMER, 1, TYPE_STRING); //name
        command_table["START_TIMER"] = START_TIMER;
        define_syntax(START_TIMER, 2, TYPE_INTEGER_VARIABLE, TYPE_DOUBLE);
        command_table["PAUSE_TIMER"] = PAUSE_TIMER;
        define_syntax(PAUSE_TIMER, 1, TYPE_INTEGER_VARIABLE);
        command_table["RESUME_TIMER"] = RESUME_TIMER;
        define_syntax(RESUME_TIMER, 1, TYPE_INTEGER_VARIABLE);
        command_table["TIMER_TO_DOUBLE"] = TIMER_TO_DOUBLE;
        define_syntax(TIMER_TO_DOUBLE, 2, TYPE_INTEGER_VARIABLE, TYPE_DOUBLE_VARIABLE);
        command_table["IF"] = IF;
        define_syntax(IF, 1, TYPE_INTEGER); //condition
        command_table["ELSE"] = ELSE;
        define_syntax(ELSE, 0); //no arguments
        command_table["ENDIF"] = ENDIF;
        define_syntax(ENDIF, 0); //no arguments
        command_table["WHILE"] = WHILE;
        define_syntax(WHILE, 1, TYPE_INTEGER); //condition
        command_table["WEND"] = WEND;
        define_syntax(WEND, 0); //no arguments
        command_table["BREAK"] = BREAK;
        define_syntax(BREAK, 0); //no arguments
        command_table["GOSUB"] = GOSUB;
        define_syntax(GOSUB, 1, TYPE_INTEGER); //line
        command_table["RETURN"] = RETURN;
        define_syntax(RETURN, 0); //no arguments
        command_table["LABEL"] = LABEL;
        define_syntax(LABEL, 1, TYPE_INTEGER_VARIABLE); //variable to hold line number
        command_table["LOAD_SPRITESHEET"] = LOAD_SPRITESHEET;
        define_syntax(LOAD_SPRITESHEET, 6, TYPE_STRING, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //filename, width, height, origin x, origin y, ID out
        command_table["LOAD_SOUND"] = LOAD_SOUND;
        define_syntax(LOAD_SOUND, 3, TYPE_STRING, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //filename, track, ID out
        command_table["LOAD_MUSIC"] = LOAD_MUSIC;
        define_syntax(LOAD_MUSIC, 2, TYPE_STRING, TYPE_INTEGER_VARIABLE); //filename, ID out
        command_table["LOAD_RECORD"] = LOAD_RECORD;
        define_syntax(LOAD_RECORD, 5, TYPE_STRING, TYPE_INTEGER, TYPE_STRING, TYPE_INTEGER, TYPE_INTEGER_VARIABLE); //filename, entry count, default name, default score, ID out
        command_table["LOAD_CONFIG"] = LOAD_CONFIG;
        define_syntax(LOAD_CONFIG, 1, TYPE_STRING); //filename
        command_table["LOAD_CONTROLS"] = LOAD_CONTROLS;
        define_syntax(LOAD_CONTROLS, 1, TYPE_STRING); //filename
        command_table["SAVE_RECORD"] = SAVE_RECORD;
        define_syntax(SAVE_RECORD, 4, TYPE_INTEGER, TYPE_STRING, TYPE_STRING, TYPE_INTEGER); //ID in, filename, name, score
        command_table["SAVE_CONFIG"] = SAVE_CONFIG;
        define_syntax(SAVE_CONFIG, 1, TYPE_STRING); //filename
        command_table["SAVE_CONTROLS"] = SAVE_CONTROLS;
        define_syntax(SAVE_CONTROLS, 1, TYPE_STRING); //filename
        command_table["GET_RECORD"] = GET_RECORD;
        define_syntax(GET_RECORD, 4, TYPE_INTEGER, TYPE_INTEGER, TYPE_STRING_VARIABLE, TYPE_INTEGER_VARIABLE); //ID in, rank, name out, score out
        command_table["QUALIFIES"] = QUALIFIES;
        define_syntax(QUALIFIES, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER); //ID in, score, rank out (0 if not qualified)
        command_table["DEFAULT_CONFIG"] = DEFAULT_CONFIG;
        define_syntax(DEFAULT_CONFIG, 0);
        command_table["DEFAULT_CONTROLS"] = DEFAULT_CONTROLS;
        define_syntax(DEFAULT_CONTROLS, 0);
        command_table["BIND_KEY"] = BIND_KEY;
        define_syntax(BIND_KEY, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER); //key, input, negative edge
        command_table["BIND_BUTTON"] = BIND_BUTTON;
        define_syntax(BIND_BUTTON, 3, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER); //button, input, negative edge
        command_table["BIND_AXIS"] = BIND_AXIS;
        define_syntax(BIND_AXIS, 5, TYPE_INTEGER, TYPE_INTEGER, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE); //axis, input, invert, dead zone, sensitivity
        command_table["UNBIND_KEY"] = UNBIND_KEY;
        define_syntax(UNBIND_KEY, 1, TYPE_INTEGER); //key
        command_table["UNBIND_BUTTON"] = UNBIND_BUTTON;
        define_syntax(UNBIND_BUTTON, 1, TYPE_INTEGER); //button
        command_table["UNBIND_AXIS"] = UNBIND_AXIS;
        define_syntax(UNBIND_AXIS, 1, TYPE_INTEGER); //axis
        command_table["UNBIND_ALL"] = UNBIND_ALL;
        define_syntax(UNBIND_ALL, 0); //no arguments
        command_table["DISPLAY_MODE"] = DISPLAY_MODE;
        define_syntax(DISPLAY_MODE, 2, TYPE_INTEGER, TYPE_INTEGER); //fullscreen, borderless
        command_table["RESOLUTION"] = RESOLUTION;
        define_syntax(RESOLUTION, 2, TYPE_INTEGER, TYPE_INTEGER); //width, height
        command_table["FRAMERATE"] = FRAMERATE;
        define_syntax(FRAMERATE, 1, TYPE_INTEGER); //framerate
        command_table["FLASHING_LIGHTS"] = FLASHING_LIGHTS;
        define_syntax(FLASHING_LIGHTS, 1, TYPE_INTEGER); //on/off
        command_table["FLASHING_LIGHTS_ENABLED"] = FLASHING_LIGHTS_ENABLED;
        define_syntax(FLASHING_LIGHTS_ENABLED, 1, TYPE_INTEGER_VARIABLE); //bool out
        command_table["MUTE"] = MUTE;
        define_syntax(MUTE, 1, TYPE_INTEGER); //on/off
        command_table["MASTER_VOLUME"] = MASTER_VOLUME;
        define_syntax(MASTER_VOLUME, 1, TYPE_INTEGER); //volume
        command_table["MUSIC_VOLUME"] = MUSIC_VOLUME;
        define_syntax(MUSIC_VOLUME, 1, TYPE_INTEGER); //volume
        command_table["SOUND_VOLUME"] = SOUND_VOLUME;
        define_syntax(SOUND_VOLUME, 1, TYPE_INTEGER); //volume
        command_table["SCALE_UI"] = SCALE_UI;
        define_syntax(SCALE_UI, 1, TYPE_DOUBLE); //scale
        command_table["SCALE_HUD"] = SCALE_HUD;
        define_syntax(SCALE_HUD, 1, TYPE_DOUBLE); //scale
        command_table["PRINT"] = PRINT;
        define_syntax(PRINT, 1, TYPE_STRING); //message
        command_table["FLIPPER"] = FLIPPER;
        define_syntax(FLIPPER, 9, TYPE_INTEGER, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_DOUBLE, TYPE_INTEGER_VARIABLE); //layer, shaft x, shaft y, angle, angular range, length, large end radius, small end radius, ID out
        //Events
        event_table["LOAD"] = LOAD;
        event_table["EVERY_TICK"] = EVERY_TICK;
        event_table["EVERY_FRAME"] = EVERY_FRAME;
        event_table["END"] = END;
        event_table["FUNCTION"] = FUNCTION;
        event_table["INPUT_DOWN"] = INPUT_DOWN;
        event_table["INPUT_UP"] = INPUT_UP;
        event_table["COLLISION"] = COLLISION;
        event_table["SWITCH_DOWN"] = SWITCH_DOWN;
        event_table["SWITCH_UP"] = SWITCH_UP;
        event_table["TILT_BOB_CONTACT"] = TILT_BOB_CONTACT;
        event_table["SPIN_CW"] = SPIN_CW;
        event_table["SPIN_CCW"] = SPIN_CCW;
        event_table["TIMER_DONE"] = TIMER_DONE;
        event_table["REEL_RESET_DONE"] = REEL_RESET_DONE;
        event_table["REEL_ROLLOVER"] = REEL_ROLLOVER;

        enum_table["INPUT_START"] = INPUT_START;
        enum_table["INPUT_COIN"] = INPUT_COIN;
        enum_table["INPUT_LEFT_FLIPPER"] = INPUT_LEFT_FLIPPER;
        enum_table["INPUT_RIGHT_FLIPPER"] = INPUT_RIGHT_FLIPPER;
        enum_table["INPUT_LEFT_UPPER_FLIPPER"] = INPUT_LEFT_UPPER_FLIPPER;
        enum_table["INPUT_RIGHT_UPPER_FLIPPER"] = INPUT_RIGHT_UPPER_FLIPPER;
        enum_table["INPUT_LEFT_SECONDARY"] = INPUT_LEFT_SECONDARY;
        enum_table["INPUT_RIGHT_SECONDARY"] = INPUT_RIGHT_SECONDARY;
        enum_table["INPUT_BALL_ELEVATOR"] = INPUT_BALL_ELEVATOR;
        enum_table["INPUT_PLUNGER"] = INPUT_PLUNGER;
        enum_table["INPUT_NUDGE"] = INPUT_NUDGE;
        enum_table["INPUT_MENU_UP"] = INPUT_MENU_UP;
        enum_table["INPUT_MENU_DOWN"] = INPUT_MENU_DOWN;
        enum_table["INPUT_MENU_LEFT"] = INPUT_MENU_LEFT;
        enum_table["INPUT_MENU_RIGHT"] = INPUT_MENU_RIGHT;
        enum_table["INPUT_MENU_CONFIRM"] = INPUT_MENU_CONFIRM;
        enum_table["INPUT_MENU_CANCEL"] = INPUT_MENU_CANCEL;
        enum_table["INPUT_CAMERA_ZOOM_IN"] = INPUT_CAMERA_ZOOM_IN;
        enum_table["INPUT_CAMERA_ZOOM_OUT"] = INPUT_CAMERA_ZOOM_OUT;
        enum_table["INPUT_CAMERA_HOME"] = INPUT_CAMERA_HOME;
        enum_table["ANALOG_PLUNGER"] = ANALOG_PLUNGER;
        enum_table["ANALOG_NUDGE_X"] = ANALOG_NUDGE_X;
        enum_table["ANALOG_NUDGE_Y"] = ANALOG_NUDGE_Y;
        enum_table["ANALOG_CAMERA_X"] = ANALOG_CAMERA_X;
        enum_table["ANALOG_CAMERA_Y"] = ANALOG_CAMERA_Y;
        enum_table["ANALOG_CAMERA_ZOOM"] = ANALOG_CAMERA_ZOOM;
        enum_table["TRACK_PHYSICAL"] = TRACK_PHYSICAL;
        enum_table["TRACK_TABLE"] = TRACK_TABLE;
        enum_table["TRACK_MENU"] = TRACK_MENU;
        enum_table["TRACK_MUSIC"] = TRACK_MUSIC;
        enum_table["ALIGN_BOTTOM_LEFT"] = ALIGN_BOTTOM_LEFT;
        enum_table["ALIGN_BOTTOM"] = ALIGN_BOTTOM;
        enum_table["ALIGN_BOTTOM_RIGHT"] = ALIGN_BOTTOM_RIGHT;
        enum_table["ALIGN_LEFT"] = ALIGN_LEFT;
        enum_table["ALIGN_CENTER"] = ALIGN_CENTER;
        enum_table["ALIGN_RIGHT"] = ALIGN_RIGHT;
        enum_table["ALIGN_TOP_LEFT"] = ALIGN_TOP_LEFT;
        enum_table["ALIGN_TOP"] = ALIGN_TOP;
        enum_table["ALIGN_TOP_RIGHT"] = ALIGN_TOP_RIGHT;
        enum_table["GRADIENT_CONSTANT"] = GRADIENT_CONSTANT;
        enum_table["GRADIENT_POLAR"] = GRADIENT_POLAR;
        enum_table["FALSE"] = 0;
        enum_table["TRUE"] = 1;
        enum_table["DISPLAY_WINDOWED"] = 0;
        enum_table["DISPLAY_FULLSCREEN"] = 1;
        enum_table["DISPLAY_BORDERED"] = 0;
        enum_table["DISPLAY_BORDERLESS"] = 1;
        enum_table["BELOW_BALL"] = 0;
        enum_table["ABOVE_BALL"] = 1;
        enum_table["STANDARD_NUDGE"] = 0;
        enum_table["UNIFORM_NUDGE"] = 1;
        //TODO (Phase VIII) key/button enums?
        //TODO gate enums?
    }
    ~Syntax()
    {
        for (auto i = arg_table.begin(); i != arg_table.end(); ++i)
        {
            delete i->second;
        }
    }
    void define_syntax(Command command, int arg_count, ...)
    {
        std::vector<DataType>* arg_list = new std::vector<DataType>;
        std::va_list passed;
        va_start(passed, arg_count);
        for (int i = 0; i < arg_count; ++i)
        {
            arg_list->push_back(static_cast<DataType>(va_arg(passed, int)));
        }
        va_end(passed);
        arg_table[command] = arg_list;
    }
    DataType get_arg_type(Command command, unsigned int arg)
    {
        if (arg < arg_table[command]->size())
        {
            return (*arg_table[command])[arg];
        } else
        {
            return TYPE_NULL;
        }
    }
};

class InputHandler
{
public:
    std::unordered_map<SDL_Keycode, std::vector<std::tuple<Input, bool, double>>> keyboard_binds; //maps keyboard inputs to game inputs with bool for negative edge and double for analog value (set when applicable)
    std::unordered_map<SDL_GameControllerButton, std::vector<std::tuple<Input, bool, double>>> gamepad_button_binds; //same binding scheme for gamepad digital inputs
    std::unordered_map<SDL_GameControllerAxis, std::vector<std::tuple<Input, bool, double, double>>> gamepad_axis_binds; //maps gamepad analog axes to inputs with bool for invert, double for deadzone, and double for sensitivity (as multiplier)
    std::pair<bool, double> input_state[static_cast<int>(INPUT_COUNT)];
    std::deque<std::pair<Event, int>> trigger_queue;
    InputHandler()
    {
        for (int i = 0; i < static_cast<int>(INPUT_COUNT) - 1; ++i)
        {
            input_state[i] = std::pair<bool, double>(false, 0.0);
        }
        default_controls();
    }
    InputHandler(char filename[])
    {
        for (int i = 0; i < static_cast<int>(INPUT_COUNT) - 1; ++i)
        {
            input_state[i] = std::pair<bool, double>(false, 0.0);
        }
        load_controls(filename);
    }
    void bind_key(SDL_Keycode key, Input input, bool negative_edge = false, double analog_value = 1.0)
    {
        keyboard_binds[key].push_back(std::tuple<Input, bool, double>(input, negative_edge, analog_value));
    }
    void bind_button(SDL_GameControllerButton button, Input input, bool negative_edge = false, double analog_value = 1.0)
    {
        gamepad_button_binds[button].push_back(std::tuple<Input, bool, double>(input, negative_edge, analog_value));
    }
    void bind_axis(SDL_GameControllerAxis axis, Input input, bool invert = false, double dead_zone = 0.0, double sensitivity = 1.0)
    {
        gamepad_axis_binds[axis].push_back(std::tuple<Input, bool, double, double>(input, invert, dead_zone, sensitivity));
    }
    void default_controls()
    {
        //TODO (Phase X) ensure portability
        keyboard_binds.clear();
        bind_key(SDLK_f, INPUT_LEFT_FLIPPER);
        bind_key(SDLK_f, INPUT_LEFT_UPPER_FLIPPER);
        bind_key(SDLK_j, INPUT_RIGHT_FLIPPER);
        bind_key(SDLK_j, INPUT_RIGHT_UPPER_FLIPPER);
        bind_key(SDLK_g, INPUT_LEFT_SECONDARY);
        bind_key(SDLK_h, INPUT_RIGHT_SECONDARY);
        bind_key(SDLK_c, INPUT_START);
        bind_key(SDLK_v, INPUT_COIN);
        bind_key(SDLK_d, ANALOG_NUDGE_X, false, 2.5);
        bind_key(SDLK_d, INPUT_NUDGE);
        bind_key(SDLK_k, ANALOG_NUDGE_X, false, -2.5);
        bind_key(SDLK_k, INPUT_NUDGE);
        bind_key(SDLK_SPACE, ANALOG_NUDGE_Y, false, -2.5);
        bind_key(SDLK_SPACE, INPUT_NUDGE);
        bind_key(SDLK_m, ANALOG_PLUNGER, false, 4.0);
        bind_key(SDLK_m, INPUT_PLUNGER, true);
        bind_key(SDLK_UP, INPUT_MENU_UP);
        bind_key(SDLK_DOWN, INPUT_MENU_DOWN);
        bind_key(SDLK_LEFT, INPUT_MENU_LEFT);
        bind_key(SDLK_RIGHT, INPUT_MENU_RIGHT);
        bind_key(SDLK_RETURN, INPUT_MENU_CONFIRM);
        bind_key(SDLK_ESCAPE, INPUT_MENU_CANCEL);
        bind_key(SDLK_KP_8, ANALOG_CAMERA_Y, false, -1.0);
        bind_key(SDLK_KP_2, ANALOG_CAMERA_Y, false, 1.0);
        bind_key(SDLK_KP_4, ANALOG_CAMERA_X, false, -1.0);
        bind_key(SDLK_KP_6, ANALOG_CAMERA_X, false, 1.0);
        bind_key(SDLK_EQUALS, INPUT_CAMERA_ZOOM_IN);
        bind_key(SDLK_KP_PLUS, INPUT_CAMERA_ZOOM_IN);
        bind_key(SDLK_MINUS, INPUT_CAMERA_ZOOM_OUT);
        bind_key(SDLK_KP_MINUS, INPUT_CAMERA_ZOOM_OUT);
        bind_key(SDLK_KP_5, INPUT_CAMERA_HOME);
        //TODO (Phase VIII) gamepad camera binds
        gamepad_button_binds.clear();
        bind_button(static_cast<SDL_GameControllerButton>(4), INPUT_LEFT_SECONDARY); //Left shoulder (SDL enum doesn't work)
        bind_button(static_cast<SDL_GameControllerButton>(5), INPUT_RIGHT_SECONDARY); //Right shoulder (SDL enum doesn't work)
        bind_button(SDL_CONTROLLER_BUTTON_START, INPUT_START);
        bind_button(SDL_CONTROLLER_BUTTON_BACK, INPUT_COIN);
        bind_button(SDL_CONTROLLER_BUTTON_DPAD_UP, INPUT_MENU_UP);
        bind_button(SDL_CONTROLLER_BUTTON_DPAD_DOWN, INPUT_MENU_DOWN);
        bind_button(SDL_CONTROLLER_BUTTON_DPAD_LEFT, INPUT_MENU_LEFT);
        bind_button(SDL_CONTROLLER_BUTTON_DPAD_RIGHT, INPUT_MENU_RIGHT);
        bind_button(SDL_CONTROLLER_BUTTON_A, INPUT_MENU_CONFIRM);
        bind_button(SDL_CONTROLLER_BUTTON_B, INPUT_MENU_CANCEL);
        bind_button(SDL_CONTROLLER_BUTTON_X, INPUT_NUDGE);
        bind_button(SDL_CONTROLLER_BUTTON_Y, INPUT_PLUNGER);
        gamepad_axis_binds.clear();
        bind_axis(SDL_CONTROLLER_AXIS_TRIGGERLEFT, INPUT_LEFT_FLIPPER, false, 0.2);
        bind_axis(SDL_CONTROLLER_AXIS_TRIGGERLEFT, INPUT_LEFT_UPPER_FLIPPER, false, 0.8);
        bind_axis(SDL_CONTROLLER_AXIS_TRIGGERRIGHT, INPUT_RIGHT_FLIPPER, false, 0.2);
        bind_axis(SDL_CONTROLLER_AXIS_TRIGGERRIGHT, INPUT_RIGHT_UPPER_FLIPPER, false, 0.8);
        bind_axis(SDL_CONTROLLER_AXIS_LEFTX, ANALOG_NUDGE_X, false, 0.0, 4.0);
        bind_axis(SDL_CONTROLLER_AXIS_LEFTY, ANALOG_NUDGE_Y, false, 0.0, 4.0);
        bind_axis(SDL_CONTROLLER_AXIS_RIGHTY, ANALOG_PLUNGER, false, 0.4, 4.0); //TODO (Phase X) test analog plunger on all tables
    }
    void load_controls(char filename[])
    {
        keyboard_binds.clear();
        gamepad_button_binds.clear();
        gamepad_axis_binds.clear();
        std::ifstream file(filename);
        std::string line;
        int mode = 0; //0 = key, 1 = button, 2 = axis
        int code = 0;
        while (std::getline(file, line))
        {
            std::string word;
            std::istringstream linestream(line);
            linestream >> word;
            if (word == "Key")
            {
                mode = 0;
                linestream >> word;
                code = std::stoi(word);
            } else if (word == "Button")
            {
                mode = 1;
                linestream >> word;
                code = std::stoi(word);
            } else if (word == "Axis")
            {
                mode = 2;
                linestream >> word;
                code = std::stoi(word);
            } else
            {
                switch (mode)
                {
                case 0:
                    {
                        Input input = static_cast<Input>(std::stoi(word));
                        linestream >> word;
                        bool negative_edge = std::stoi(word);
                        linestream >> word;
                        double analog_value = std::stod(word);
                        bind_key(code, input, negative_edge, analog_value);
                    }
                    break;
                case 1:
                    {
                        Input input = static_cast<Input>(std::stoi(word));
                        linestream >> word;
                        bool negative_edge = std::stoi(word);
                        linestream >> word;
                        double analog_value = std::stod(word);
                        bind_button(static_cast<SDL_GameControllerButton>(code), input, negative_edge, analog_value);
                    }
                    break;
                case 2:
                    {
                        Input input = static_cast<Input>(std::stoi(word));
                        linestream >> word;
                        bool invert = std::stoi(word);
                        linestream >> word;
                        double dead_zone = std::stod(word);
                        linestream >> word;
                        double sensitivity = std::stod(word);
                        bind_axis(static_cast<SDL_GameControllerAxis>(code), input, invert, dead_zone, sensitivity);
                    }
                    break;
                }
            }
        }
    }
    void save_controls(char filename[])
    {
        std::ofstream file(filename);
        for (auto key_vector = keyboard_binds.begin(); key_vector != keyboard_binds.end(); ++key_vector)
        {
            file << "Key " << key_vector->first << std::endl;
            for (auto bind_tuple = key_vector->second.begin(); bind_tuple != key_vector->second.end(); ++bind_tuple)
            {
                file << std::to_string(static_cast<int>(std::get<0>(*bind_tuple))) << " " << std::to_string(std::get<1>(*bind_tuple)) << " " << std::to_string(std::get<2>(*bind_tuple)) << std::endl;
            }
        }
        for (auto button_vector = gamepad_button_binds.begin(); button_vector != gamepad_button_binds.end(); ++button_vector)
        {
            file << "Button " << button_vector->first << std::endl;
            for (auto bind_tuple = button_vector->second.begin(); bind_tuple != button_vector->second.end(); ++bind_tuple)
            {
                file << std::to_string(static_cast<int>(std::get<0>(*bind_tuple))) << " " << std::to_string(std::get<1>(*bind_tuple)) << " " << std::to_string(std::get<2>(*bind_tuple)) << std::endl;
            }
        }
        for (auto axis_vector = gamepad_axis_binds.begin(); axis_vector != gamepad_axis_binds.end(); ++axis_vector)
        {
            file << "Axis " << axis_vector->first << std::endl;
            for (auto bind_tuple = axis_vector->second.begin(); bind_tuple != axis_vector->second.end(); ++bind_tuple)
            {
                file << std::to_string(static_cast<int>(std::get<0>(*bind_tuple))) << " " << std::to_string(std::get<1>(*bind_tuple)) << " " << std::to_string(std::get<2>(*bind_tuple)) << " " << std::to_string(std::get<3>(*bind_tuple)) << std::endl;
            }
        }
    }
    void process_event(SDL_Event* e)
    {
        switch(e->type)
        {
        case SDL_KEYDOWN:
            {
                if (keyboard_binds.find(e->key.keysym.sym) != keyboard_binds.end() && !e->key.repeat)
                {
                    auto bind_list = keyboard_binds[e->key.keysym.sym];
                    for (unsigned int i = 0; i < bind_list.size(); ++i)
                    {
                        auto bind_data = bind_list[i];
                        Input input = std::get<0>(bind_data);
                        bool negative_edge = std::get<1>(bind_data);
                        if (!negative_edge)
                        {
                            input_state[input].first = true;
                            input_state[input].second = std::get<2>(bind_data);
                            trigger_queue.push_back(std::pair<Event, int>(INPUT_DOWN, static_cast<int>(input)));
                        } else
                        {
                            input_state[input].first = false;
                            input_state[input].second = 0.0;
                            trigger_queue.push_back(std::pair<Event, int>(INPUT_UP, static_cast<int>(input)));
                        }
                    }
                }
            }
            break;
        case SDL_KEYUP:
            {
                if (keyboard_binds.find(e->key.keysym.sym) != keyboard_binds.end() && !e->key.repeat)
                {
                    auto bind_list = keyboard_binds[e->key.keysym.sym];
                    for (unsigned int i = 0; i < bind_list.size(); ++i)
                    {
                        auto bind_data = bind_list[i];
                        Input input = std::get<0>(bind_data);
                        bool negative_edge = std::get<1>(bind_data);
                        if (negative_edge)
                        {
                            input_state[input].first = true;
                            input_state[input].second = std::get<2>(bind_data);
                            trigger_queue.push_back(std::pair<Event, int>(INPUT_DOWN, static_cast<int>(input)));
                        } else
                        {
                            input_state[input].first = false;
                            input_state[input].second = 0.0;
                            trigger_queue.push_back(std::pair<Event, int>(INPUT_UP, static_cast<int>(input)));
                        }
                    }
                }
            }
            break;
        case SDL_JOYBUTTONDOWN:
            {
                if (gamepad_button_binds.find(static_cast<SDL_GameControllerButton>(e->cbutton.button)) != gamepad_button_binds.end())
                {
                    auto bind_list = gamepad_button_binds[static_cast<SDL_GameControllerButton>(e->cbutton.button)];
                    for (unsigned int i = 0; i < bind_list.size(); ++i)
                    {
                        auto bind_data = bind_list[i];
                        Input input = std::get<0>(bind_data);
                        bool negative_edge = std::get<1>(bind_data);
                        if (!negative_edge)
                        {
                            input_state[input].first = true;
                            input_state[input].second = std::get<2>(bind_data);
                            trigger_queue.push_back(std::pair<Event, int>(INPUT_DOWN, static_cast<int>(input)));
                        } else
                        {
                            input_state[input].first = false;
                            input_state[input].second = 0.0;
                            trigger_queue.push_back(std::pair<Event, int>(INPUT_UP, static_cast<int>(input)));
                        }
                    }
                }
            }
            break;
        case SDL_JOYBUTTONUP:
            {
                if (gamepad_button_binds.find(static_cast<SDL_GameControllerButton>(e->cbutton.button)) != gamepad_button_binds.end())
                {
                    auto bind_list = gamepad_button_binds[static_cast<SDL_GameControllerButton>(e->cbutton.button)];
                    for (unsigned int i = 0; i < bind_list.size(); ++i)
                    {
                        auto bind_data = bind_list[i];
                        Input input = std::get<0>(bind_data);
                        bool negative_edge = std::get<1>(bind_data);
                        if (negative_edge)
                        {
                            input_state[input].first = true;
                            input_state[input].second = std::get<2>(bind_data);
                            trigger_queue.push_back(std::pair<Event, int>(INPUT_DOWN, static_cast<int>(input)));
                        } else
                        {
                            input_state[input].first = false;
                            input_state[input].second = 0.0;
                            trigger_queue.push_back(std::pair<Event, int>(INPUT_UP, static_cast<int>(input)));
                        }
                    }
                }
            }
            break;
        case SDL_JOYAXISMOTION:
            {
                if (gamepad_axis_binds.find(static_cast<SDL_GameControllerAxis>(e->caxis.axis)) != gamepad_axis_binds.end())
                {
                    auto bind_list = gamepad_axis_binds[static_cast<SDL_GameControllerAxis>(e->caxis.axis)];
                    for (unsigned int i = 0; i < bind_list.size(); ++i)
                    {
                        auto bind_data = bind_list[i];
                        Input input = std::get<0>(bind_data);
                        bool invert = std::get<1>(bind_data);
                        double dead_zone = std::get<2>(bind_data);
                        double sensitivity = std::get<3>(bind_data);
                        double value = e->caxis.value / 32768.0 * (invert ? -1.0 : 1.0);
                        if (e->caxis.axis == SDL_CONTROLLER_AXIS_TRIGGERLEFT || e->caxis.axis == SDL_CONTROLLER_AXIS_TRIGGERRIGHT)
                        {
                            value = (value + 1.0) / 2.0;
                        }
                        value *= sensitivity;
                        if (std::abs(value) < dead_zone - ANALOG_INPUT_EPSILON)
                        {
                            if (input_state[input].first)
                            {
                                trigger_queue.push_back(std::pair<Event, int>(INPUT_UP, static_cast<int>(input)));
                            }
                            input_state[input].first = false;
                            input_state[input].second = 0.0;
                        } else if (std::abs(value) > dead_zone)
                        {
                            if (!input_state[input].first)
                            {
                                trigger_queue.push_back(std::pair<Event, int>(INPUT_DOWN, static_cast<int>(input)));
                            }
                            input_state[input].first = true;
                            input_state[input].second = value;
                        }
                    }
                }
            }
            break;
        }
    }
};

class Record
{
public:
    std::vector<std::pair<std::string, int>> scores;
    Record(char filename[], int entry_count, std::string default_name, int default_score)
    {
        std::ifstream file(filename);
        if (file.good())
        {
            load(filename);
        } else
        {
            for (int i = 0; i < entry_count; ++i)
            {
                scores.push_back(std::pair<std::string, int>(default_name, default_score));
            }
        }
    }
    void load(char filename[])
    {
        scores.clear();
        std::ifstream file(filename);
        std::string line;
        while (std::getline(file, line))
        {
            std::string word;
            std::istringstream linestream(line);
            linestream >> word;
            std::string name = word;
            linestream >> word;
            int score = std::stoi(word);
            scores.push_back(std::pair<std::string, int>(name, score));
        }
    }
    void save(char filename[])
    {
        std::ofstream file(filename);
        for (unsigned int i = 0; i < scores.size(); ++i)
        {
            file << scores[i].first << " " << scores[i].second << std::endl;
        }
    }
    int qualifies(int score)
    {
        if (score <= scores.back().second)
        {
            return 0;
        }
        for (unsigned int i = 0; i < scores.size(); i--)
        {
            if (score > scores[i].second)
            {
                return i + 1;
            }
        }
        return 0;
    }
    void add_record(std::string name, int score)
    {
        /*
        if (!qualifies(score))
        {
            return;
        }
        for (unsigned int i = 0; i < scores.size(); ++i)
        {
            if (score > scores[i].second)
            {
                scores.insert(scores.begin() + i, std::pair<std::string, int>(name, score));
                scores.pop_back();
                break;
            }
        }
        */
        //TODO (Phase X) test binary search (including tied scores)
        int high = scores.size();
        int low = 0;
        while (high > low)
        {
            int middle = (high + low) / 2;
            int a = scores[middle].second;
            int b = score;
            if (a >= b)
            {
                low = middle;
            } else
            {
                high = middle;
            }
        }
        if (static_cast<unsigned int>(low) < scores.size())
        {
            scores.insert(scores.begin() + low, std::pair<std::string, int>(name, score));
            scores.pop_back();
        }
    }
};

class Instruction
{
public:
    Command command;
    std::vector<std::tuple<DataType, int, int>> variable_table; //tuple is data type, variable index, argument index
    std::vector<int> int_args; //Note: stores integers and indices of variables for writing out
    std::vector<double> double_args;
    std::vector<std::string> string_args;
};

class Layer
{
public:
    Point top_left = Point(0.0, 0.0);
    Point dimensions = Point(20.25, 42.0);
    int next_layer = -1;
    std::vector<Collider*> colliders;
    std::vector<MobileCollider*> mobile_colliders;
    std::deque<std::pair<Ball*, SpriteInstance*>> balls;
    std::vector<Collider*>** partitioning_grid = nullptr;
    std::vector<Magnet*> magnets;
    std::vector<Magnet*>** magnet_partitioning_grid = nullptr; //TODO (Phase X) test partitioning grid system against errors
    int partitioning_grid_width = 4;
    int partitioning_grid_height = 4;
    GradientMode gradient_mode = GRADIENT_CONSTANT;
    Point gradient_vector = Point(0.0, 40.136165); //adjust depending on slope of table
    //TODO add vector offset applied after gradient vector is calculated
    Point center = Point(dimensions.x / 2.0, dimensions.y / 2.0);
    Point ball_return = Point(center.x, center.y);
    std::vector<SpriteInstance*> sprites_below_ball;
    std::vector<SpriteInstance*> sprites_above_ball;
    std::vector<Backglass*> backglasses;
    std::deque<std::pair<Event, int>> trigger_queue;
    std::unordered_set<int> switches_down;
    double spin_friction = 2.5;
    double rolling_friction = 2.75;
    int id = -1;
    double highest_radius = 0.0;
    Layer()
    {
        //TODO (Phase X) test various widths and heights (including memory leak testing on unload)
        partitioning_grid = new std::vector<Collider*>*[partitioning_grid_width];
        magnet_partitioning_grid = new std::vector<Magnet*>*[partitioning_grid_width];
        for (int column = 0; column < partitioning_grid_width; ++column)
        {
            partitioning_grid[column] = new std::vector<Collider*>[partitioning_grid_height];
            magnet_partitioning_grid[column] = new std::vector<Magnet*>[partitioning_grid_height];
        }
    }
    ~Layer()
    {
        //TODO (Phase X) check if partitioning grid sub-arrays must be deleted before main array to avoid memory leak
        delete[] partitioning_grid;
        delete[] magnet_partitioning_grid;
        for (unsigned int i = 0; i < balls.size(); ++i)
        {
            delete balls[i].second;
        }
        for (unsigned int i = 0; i < sprites_above_ball.size(); ++i)
        {
            delete sprites_above_ball[i];
        }
        for (unsigned int i = 0; i < sprites_below_ball.size(); ++i)
        {
            delete sprites_below_ball[i];
        }
    }
    int partitioning_grid_row(double y)
    {
        return std::max(std::min(static_cast<int>(y * partitioning_grid_height / dimensions.y), partitioning_grid_height - 1), 0);
    }
    int partitioning_grid_column(double x)
    {
        return std::max(std::min(static_cast<int>(x * partitioning_grid_width / dimensions.x), partitioning_grid_width - 1), 0);
    }
    void add_ball(Ball* b, Sprite* s = nullptr)
    {
        if (s)
        {
            if (b->sprite_rotation)
            {
                balls.push_back(std::pair<Ball*, SpriteInstance*>(b, new SpriteInstance(s, b->position.x, b->position.y, b->angular_displacement)));
            } else
            {
                balls.push_back(std::pair<Ball*, SpriteInstance*>(b, new SpriteInstance(s, b->position.x, b->position.y, 0.0)));
            }
        } else
        {
            balls.push_back(std::pair<Ball*, SpriteInstance*>(b, nullptr));
        }
        if (b->radius > highest_radius)
        {
            highest_radius = b->radius;
        }
    }
    void pass_ball(Ball* ball, Layer* destination)
    {
        int index = -1;
        for (unsigned int i = 0; i < balls.size(); ++i)
        {
            if (ball == balls[i].first)
            {
                index = i;
                break;
            }
        }
        if (index != -1 && destination && !(destination == this))
        {
            if (balls[index].second)
            {
                destination->add_ball(balls[index].first, balls[index].second->sprite);
            } else
            {
                destination->add_ball(balls[index].first, nullptr);
            }
            balls[index].first->layer = destination->id;
            balls.erase(balls.begin() + index);
        }
    }
    SpriteInstance* add_sprite(Sprite* s, bool above_ball = false, MobileCollider* mc = nullptr)
    {
        SpriteInstance* instance = new SpriteInstance(s);
        if (mc)
        {
            instance->sync(mc);
            instance->synced_collider = mc;
        }
        if (above_ball)
        {
            sprites_above_ball.push_back(instance);
        } else
        {
            sprites_below_ball.push_back(instance);
        }
        return instance;
    }
    SpriteInstance* add_sprite(Sprite* s, double x, double y, bool above_ball = false, double angle = 0.0)
    {
        SpriteInstance* instance = new SpriteInstance(s, x, y, angle);
        if (above_ball)
        {
            sprites_above_ball.push_back(instance);
        } else
        {
            sprites_below_ball.push_back(instance);
        }
        return instance;
    }
    void add_collider(Collider* c)
    {
        if (c->id == -1)
        {
            std::cout << "WARNING: Collider ID not set" << std::endl;
        }
        int first_row = partitioning_grid_row(c->bounding_box_top_left.y);
        int last_row = partitioning_grid_row(c->bounding_box_bottom_right.y);
        int first_column = partitioning_grid_column(c->bounding_box_top_left.x);
        int last_column = partitioning_grid_column(c->bounding_box_bottom_right.x);
        for (int row = first_row; row <= last_row; ++row)
        {
            for (int column = first_column; column <= last_column; ++column)
            {
                partitioning_grid[row][column].push_back(c);
            }
        }
        colliders.push_back(c);
        if (c->mobile)
        {
            mobile_colliders.push_back(dynamic_cast<MobileCollider*>(c));
        }
    }
    void add_magnet(Magnet* m)
    {
        int first_row = partitioning_grid_row(m->position.y - m->radius);
        int last_row = partitioning_grid_row(m->position.y + m->radius);
        int first_column = partitioning_grid_column(m->position.x - m->radius);
        int last_column = partitioning_grid_column(m->position.x + m->radius);
        for (int row = first_row; row <= last_row; ++row)
        {
            for (int column = first_column; column <= last_column; ++column)
            {
                magnet_partitioning_grid[row][column].push_back(m);
            }
        }
        magnets.push_back(m);
    }
    Point gradient_at_point(Point* p)
    {
        switch (gradient_mode)
        {
        case GRADIENT_CONSTANT:
            return gradient_vector;
            break;
        case GRADIENT_POLAR:
            {
                Point center_to_point = Point(p->x - center.x, p->y - center.y);
                if (center_to_point.square_magnitude() == 0.0)
                {
                    return Point(0.0, 0.0);
                }
                double center_to_point_angle = std::atan2(center_to_point.y, center_to_point.x);
                center_to_point_angle += PI / 2.0; //Corrects for standard (down = towards center)
                Point rotated_gradient_vector = Point(gradient_vector.x * std::cos(center_to_point_angle) - gradient_vector.y * std::sin(center_to_point_angle), gradient_vector.x * std::sin(center_to_point_angle) + gradient_vector.y * std::cos(center_to_point_angle));
                return rotated_gradient_vector;
            }
            break;
        default:
            return Point(0.0, 0.0);
            break;
        }
    }
    void simulate()
    {
        std::unordered_set<int> switches_down_this_tick;
        for (unsigned int i = 0; i < mobile_colliders.size(); ++i)
        {
            mobile_colliders[i]->simulate();
        }
        //Note: benchmarking indicates insertion sort is marginally faster than std::sort
        for (unsigned int i = 1; i < balls.size(); ++i)
        {
            if (balls[i].first->position.y > balls[i - 1].first->position.y)
            {
                for (int j = i - 1; j >= 0; --j)
                {
                    if (balls[j + 1].first->position.y > balls[j].first->position.y)
                    {
                        std::swap(balls[j], balls[j + 1]);
                    } else
                    {
                        break;
                    }
                }
            }
        }
        for (unsigned int i = 0; i + 1 < balls.size(); ++i)
        {
            Ball* ball_a = balls[i].first;
            for (unsigned int j = i + 1; j < balls.size(); ++j)
            {
                Ball* ball_b = balls[j].first;
                if (ball_b->position.y - ball_a->position.y > 2.0 * highest_radius)
                {
                    break;
                }
                ball_a->is_touching(ball_b, true, 2.0); //TODO (Phase X) test with different values (1.5 may be realistic), set bounce coefficient by script?
            }
        }
        for (unsigned int i = 0; i < balls.size(); ++i)
        {
            Ball* b = balls[i].first;
            b->friction_multiplier = 1.0;
            int first_row = partitioning_grid_row(b->position.y - b->radius);
            int last_row = partitioning_grid_row(b->position.y + b->radius);
            int first_column = partitioning_grid_column(b->position.x - b->radius);
            int last_column = partitioning_grid_column(b->position.x + b->radius);
            for (int row = first_row; row <= last_row; ++row)
            {
                for (int column = first_column; column <= last_column; ++column)
                {
                    std::vector<Magnet*>* cell_magnets = &(magnet_partitioning_grid[row][column]);
                    for (unsigned int j = 0; j < cell_magnets->size(); ++j)
                    {
                        (*cell_magnets)[j]->attract(b);
                    }
                }
            }
            Point acceleration_vector = gradient_at_point(&b->position);
            if (b->friction_multiplier < 0.0)
            {
                b->friction_multiplier = 0.0;
            }
            double effective_rolling_friction = rolling_friction * b->friction_multiplier;
            if (b->velocity.magnitude() * SIM_TIME_PER_PHYSICS_TICK > effective_rolling_friction * SIM_TIME_PER_PHYSICS_TICK * SIM_TIME_PER_PHYSICS_TICK)
            {
                Point velocity_unit = Point(b->velocity.x / b->velocity.magnitude(), b->velocity.y / b->velocity.magnitude());
                acceleration_vector.x -= effective_rolling_friction * velocity_unit.x;
                acceleration_vector.y -= effective_rolling_friction * velocity_unit.y;
            } else
            {
                acceleration_vector.x -= b->velocity.x / SIM_TIME_PER_PHYSICS_TICK;
                acceleration_vector.y -= b->velocity.y / SIM_TIME_PER_PHYSICS_TICK;
            }
            if (b->velocity.magnitude() / 1000.0 * MS_PER_FRAME / PHYSICS_TICKS_PER_FRAME > b->radius / 4.0)
            {
                std::cout << "WARNING: fast ball (" << b->velocity.magnitude() << " inches per second)" << std::endl; //TODO (Phase X) remove
            }
            b->velocity.x += acceleration_vector.x * SIM_TIME_PER_PHYSICS_TICK;
            b->velocity.y += acceleration_vector.y * SIM_TIME_PER_PHYSICS_TICK;
            b->position.x += b->velocity.x * SIM_TIME_PER_PHYSICS_TICK;
            b->position.y += b->velocity.y * SIM_TIME_PER_PHYSICS_TICK;
            double energy = b->velocity.square_magnitude();
            if (b->roll_sound)
            {
                int volume = energy * 128 / b->max_roll_volume_speed / b->max_roll_volume_speed;
                if (volume > 128)
                {
                    volume = 128;
                }
                b->roll_sound->set_loop_volume(volume);
            }
            double effective_spin_friction = spin_friction * b->friction_multiplier;
            if (std::abs(b->angular_velocity) < effective_spin_friction * SIM_TIME_PER_PHYSICS_TICK)
            {
                b->angular_velocity = 0.0;
            } else if (b->angular_velocity > 0.0)
            {
                b->angular_velocity -= effective_spin_friction * SIM_TIME_PER_PHYSICS_TICK;
            } else
            {
                b->angular_velocity += effective_spin_friction * SIM_TIME_PER_PHYSICS_TICK;
            }
            b->angular_displacement += b->angular_velocity * SIM_TIME_PER_PHYSICS_TICK;
            if (b->position.x < top_left.x || b->position.x > top_left.x + dimensions.x || b->position.y < top_left.y || b->position.y > top_left.y + dimensions.y)
            {
                if (next_layer != -1)
                {
                    //TODO (Phase X) test layer switch on OOB
                    b->layer_pending = next_layer;
                } else
                {
                    b->velocity.x = 0.0;
                    b->velocity.y = 0.0;
                    b->position.x = ball_return.x;
                    b->position.y = ball_return.y;
                }
                continue;
            }
            for (int row = first_row; row <= last_row; ++row)
            {
                for (int column = first_column; column <= last_column; ++column)
                {
                    std::vector<Collider*>* cell_colliders = &(partitioning_grid[row][column]);
                    for (unsigned int j = 0; j < cell_colliders->size(); ++j)
                    {
                        Point unit_normal = Point(0.0, 0.0);
                        Collider* c = (*cell_colliders)[j];
                        if (c->is_touching(b, c->solid, &unit_normal))
                        {
                            if (c->solid)
                            {
                                double bounce_coefficient = c->bounce_coefficient;
                                LineCollider* lc = dynamic_cast<LineCollider*>(c);
                                if (lc && lc->band_bounce)
                                {
                                    Point a_to_b = Point(lc->b.x - lc->a.x, lc->b.y - lc->a.y);
                                    double a_to_b_magnitude = a_to_b.magnitude();
                                    if (a_to_b_magnitude != 0.0)
                                    {
                                        Point a_to_b_unit = Point(a_to_b.x / a_to_b_magnitude, a_to_b.y / a_to_b_magnitude);
                                        Point a_to_ball = Point(b->position.x - lc->a.x, b->position.y - lc->a.y);
                                        double ball_bias = a_to_b_unit.dot_product(a_to_ball) / a_to_b_magnitude;
                                        double bias_match = ball_bias > 0.5 ? (1.0 - ball_bias) / 0.5 : ball_bias / 0.5;
                                        double root_bias_match = std::sqrt(bias_match);
                                        bounce_coefficient = root_bias_match * lc->bounce_coefficient + (1.0 - root_bias_match) * lc->outer_bounce_coefficient;
                                    }
                                }
                                if (c->has_detection_unit_normal)
                                {
                                    b->bounce(&unit_normal, bounce_coefficient, c->friction_coefficient, c->impact_sound, &c->detection_unit_normal);
                                } else
                                {
                                    b->bounce(&unit_normal, bounce_coefficient, c->friction_coefficient, c->impact_sound);
                                }
                                //TODO (Phase X) look into possible bugs resulting from this method of collision detection
                                if (b->detection_speed >= c->sensitivity)
                                {
                                    trigger_queue.push_back(std::pair<Event, int>(COLLISION, c->id));
                                }
                            } else
                            {
                                if (switches_down.find(c->switch_id) == switches_down.end())
                                {
                                    trigger_queue.push_back(std::pair<Event, int>(SWITCH_DOWN, c->switch_id));
                                    switches_down.insert(c->switch_id);
                                }
                                if (switches_down_this_tick.find(c->switch_id) == switches_down_this_tick.end())
                                {
                                    switches_down_this_tick.insert(c->switch_id);
                                }
                            }
                            if (c->portal)
                            {
                                b->layer_pending = c->destination_layer;
                            }
                        }
                    }
                }
            }
        }
        for (auto i = switches_down.begin(); i != switches_down.end(); ++i)
        {
            if (switches_down_this_tick.find(*i) == switches_down_this_tick.end())
            {
                switches_down.erase(*i);
                trigger_queue.push_back(std::pair<Event, int>(SWITCH_UP, (*i)));
            }
        }
        for (unsigned int i = 0; i < balls.size(); ++i)
        {
            if (balls[i].second)
            {
                balls[i].second->pass_time(SIM_TIME_PER_PHYSICS_TICK * 1000.0);
            }
        }
        for (unsigned int i = 0; i < sprites_below_ball.size(); ++i)
        {
            sprites_below_ball[i]->pass_time(SIM_TIME_PER_PHYSICS_TICK * 1000.0);
        }
        for (unsigned int i = 0; i < sprites_above_ball.size(); ++i)
        {
            sprites_above_ball[i]->pass_time(SIM_TIME_PER_PHYSICS_TICK * 1000.0);
        }
        for (unsigned int i = 0; i < backglasses.size(); ++i)
        {
            //TODO stop score reel events from being framerate-dependent
            backglasses[i]->pass_time(SIM_TIME_PER_PHYSICS_TICK * 1000.0);
            pump_triggers(backglasses[i]);
        }
    }
    void pump_triggers(Backglass* backglass)
    {
        while (backglass->trigger_queue.size())
        {
            trigger_queue.push_back(backglass->trigger_queue.front());
            backglass->trigger_queue.pop_front();
        }
    }
    void add_backglass(Backglass* backglass)
    {
        backglasses.push_back(backglass);
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        if (view->sprites)
        {
            for (unsigned int i = 0; i < sprites_below_ball.size(); ++i)
            {
                if (sprites_below_ball[i]->synced_collider)
                {
                    sprites_below_ball[i]->sync(sprites_below_ball[i]->synced_collider);
                }
                sprites_below_ball[i]->render(renderer, view);
            }
        }
        for (unsigned int i = 0; i < backglasses.size(); ++i)
        {
            backglasses[i]->render(renderer, view);
        }
        for (unsigned int i = 0; i < balls.size(); ++i)
        {
            if (balls[i].second)
            {
                balls[i].second->sync(balls[i].first);
                balls[i].second->render(renderer, view);
            }
            balls[i].first->render(renderer, view);
        }
        for (unsigned int i = 0; i < sprites_above_ball.size(); ++i)
        {
            if (sprites_above_ball[i]->synced_collider)
            {
                sprites_above_ball[i]->sync(sprites_above_ball[i]->synced_collider);
            }
            sprites_above_ball[i]->render(renderer, view);
        }
        if (view->debug)
        {
            SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
            Point screen_center = Point(0.0, 0.0);
            view->to_screen_space(&center, &screen_center);
            SDL_RenderDrawPoint(renderer, screen_center.x, screen_center.y);
            SDL_Rect bounds;
            Point top_left_pixel = Point(0.0, 0.0);
            Point bottom_right = Point(top_left.x + dimensions.x, top_left.y + dimensions.y);
            Point bottom_right_pixel = Point(0.0, 0.0);
            view->to_screen_space(&top_left, &top_left_pixel);
            view->to_screen_space(&bottom_right, &bottom_right_pixel);
            bounds.x = top_left_pixel.x;
            bounds.y = top_left_pixel.y;
            bounds.w = bottom_right_pixel.x - top_left_pixel.x;
            bounds.h = bottom_right_pixel.y - top_left_pixel.y;
            SDL_SetRenderDrawColor(renderer, 128, 128, 128, 255);
            SDL_RenderDrawRect(renderer, &bounds);
            /*
            SDL_Rect cell;
            cell.w = bounds.w / partitioning_grid_width;
            cell.h = bounds.h / partitioning_grid_height;
            SDL_SetRenderDrawColor(renderer, 64, 64, 64, 255);
            for (int row = 0; row < partitioning_grid_height; ++row)
            {
                for (int column = 0; column < partitioning_grid_width; ++column)
                {
                    cell.x = bounds.x + dimensions.x * view->zoom * PIXELS_PER_INCH * column / partitioning_grid_width;
                    cell.y = bounds.y + dimensions.y * view->zoom * PIXELS_PER_INCH * row / partitioning_grid_height;
                    SDL_RenderDrawRect(renderer, &cell);
                }
            }
            //*/
            for (unsigned int i = 0; i < magnets.size(); ++i)
            {
                magnets[i]->render(renderer, view);
            }
            for (unsigned int i = 0; i < colliders.size(); ++i)
            {
                colliders[i]->render(renderer, view);
            }
        }
    }
};

class Spinner
{
public:
    double center_mass_radius = 0.5;
    double contact_point_radius = 0.5;
    double damping_coefficient = 0.50;
    double gravity = 386.0886;
    double angle = 0.0; //measured clockwise from straight down
    double angular_velocity = 0.0;
    Point unit_normal = Point(0.0, 0.0);
    Point center_mass = Point(0.0, 0.0);
    std::deque<std::pair<Event, int>> trigger_queue;
    int id = -1;
    SpriteInstance* synced_sprite = nullptr;
    int frames_per_rotation = 32;
    void simulate()
    {
        angular_velocity -= gravity * std::sin(angle) * center_mass_radius * SIM_TIME_PER_PHYSICS_TICK;
        angular_velocity *= std::pow(1.0 - damping_coefficient, SIM_TIME_PER_PHYSICS_TICK);
        angle += angular_velocity * SIM_TIME_PER_PHYSICS_TICK;
        while (angle < 0.0)
        {
            angle += 2.0 * PI;
        }
        while (angle > 2.0 * PI)
        {
            angle -= 2.0 * PI;
        }
        if (synced_sprite)
        {
            synced_sprite->frame = static_cast<int>(std::round(angle / (2 * PI / frames_per_rotation))) % frames_per_rotation;
        }
        if (angular_velocity > 0.0 && angle > PI && angle - angular_velocity * SIM_TIME_PER_PHYSICS_TICK < PI)
        {
            trigger_queue.push_back(std::pair<Event, int>(SPIN_CW, id));
        }
        if (angular_velocity < 0.0 && angle < PI && angle - angular_velocity * SIM_TIME_PER_PHYSICS_TICK > PI)
        {
            trigger_queue.push_back(std::pair<Event, int>(SPIN_CCW, id));
        }
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        if (view->debug)
        {
            Point center = Point(32.0 + 64.0 * id, 32.0);
            Point offset = Point(-std::sin(angle) * 32.0, std::cos(angle) * 32.0);
            SDL_SetRenderDrawColor(renderer, 128, 128, 128, 255);
            SDL_RenderDrawLineF(renderer, center.x, center.y, center.x + offset.x, center.y + offset.y);
        }
    }
};

void rip_spinner(Spinner* spinner, double contact_speed)
{
    spinner->angular_velocity = contact_speed / spinner->contact_point_radius;
}

class TiltBob
{
public:
    Point position = Point(0.0, 0.0);
    Point bob_position = Point(0.0, 0.0);
    Point velocity = Point(0.0, 0.0);
    double bounce_coefficient = 1.55;
    double damping_coefficient = 0.25;
    double gravity = 100.0;
    double tolerance = 0.375;
    int id = -1;
    std::deque<std::pair<Event, int>> trigger_queue;
    void simulate()
    {
        velocity.x -= bob_position.x * gravity * SIM_TIME_PER_PHYSICS_TICK;
        velocity.y -= bob_position.y * gravity * SIM_TIME_PER_PHYSICS_TICK;
        velocity.x *= std::pow(1.0 - damping_coefficient, SIM_TIME_PER_PHYSICS_TICK);
        velocity.y *= std::pow(1.0 - damping_coefficient, SIM_TIME_PER_PHYSICS_TICK);
        bob_position.x += velocity.x * SIM_TIME_PER_PHYSICS_TICK;
        bob_position.y += velocity.y * SIM_TIME_PER_PHYSICS_TICK;
        if (bob_position.square_magnitude() > tolerance * tolerance)
        {
            Point unit_bob_position = Point(bob_position.x / bob_position.magnitude(), bob_position.y / bob_position.magnitude());
            bob_position.x = unit_bob_position.x * tolerance;
            bob_position.y = unit_bob_position.y * tolerance;
            velocity.x -= velocity.x * bounce_coefficient;
            velocity.y -= velocity.y * bounce_coefficient;
            trigger_queue.push_back(std::pair<Event, int>(TILT_BOB_CONTACT, id));
        }
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        if (view->debug)
        {
            Point center = Point(0.0, 0.0);
            view->to_screen_space(&position, &center);
            int precision = 32;
            Point last_offset = Point(0.0, 0.0);
            for (int i = 0; i <= precision; ++i)
            {
                double angle = 2 * PI * i / precision;
                Point offset = Point(std::cos(angle) * tolerance * PIXELS_PER_INCH * view->zoom, std::sin(angle) * tolerance * PIXELS_PER_INCH * view->zoom);
                if (i)
                {
                    SDL_SetRenderDrawColor(renderer, 128, 128, 128, 255);
                    SDL_RenderDrawLine(renderer, center.x + last_offset.x, center.y + last_offset.y, center.x + offset.x, center.y + offset.y);
                }
                last_offset = offset;
            }
            Point bob_screen_position = Point(0.0, 0.0);
            Point bob_real_position = Point(position.x + bob_position.x, position.y + bob_position.y);
            view->to_screen_space(&bob_real_position, &bob_screen_position);
            SDL_SetRenderDrawColor(renderer, 128, 128, 128, 255);
            SDL_RenderDrawPoint(renderer, bob_screen_position.x, bob_screen_position.y);
        }
    }
};

void convert_nudge(Point* nudge_vector, Point* position, Point* nudge_origin, Point* table_center_mass, Point* positional_nudge_vector)
{
    //Note: assumes center_to_origin is parallel to y axis
    Point center_to_origin = Point(nudge_origin->x - table_center_mass->x, nudge_origin->y - table_center_mass->y);
    Point center_to_position = Point(position->x - table_center_mass->x, position->y - table_center_mass->y);
    double center_to_origin_magnitude = center_to_origin.magnitude();
    double center_to_position_magnitude = center_to_position.magnitude();
    double scale = center_to_position_magnitude / center_to_origin_magnitude;
    Point center_to_origin_unit = Point(center_to_origin.x / center_to_origin_magnitude, center_to_origin.y / center_to_origin_magnitude);
    Point center_to_position_unit = Point(center_to_position.x / center_to_position_magnitude, center_to_position.y / center_to_position_magnitude);
    double angle = std::atan2(center_to_position_unit.y, center_to_position_unit.x) - std::atan2(center_to_origin_unit.y, center_to_origin_unit.x);
    Point vector_out = Point(nudge_vector->x * scale * std::cos(angle), nudge_vector->y + nudge_vector->x * scale * std::sin(angle));
    positional_nudge_vector->x = vector_out.x;
    positional_nudge_vector->y = vector_out.y;
}

class Table
{
public:
    Point dimensions = Point(20.25, 42.0);
    std::vector<Layer*> layers;
    Backglass* backglass = nullptr;
    std::deque<std::pair<Event, int>> trigger_queue;
    std::vector<TiltBob*> tilt_bobs;
    std::vector<Spinner*> spinners;
    Point current_nudge = Point(0.0, 0.0);
    double nudge_duration = 0.0;
    double global_time = 0.0;
    bool uniform_nudge = false;
    ~Table()
    {
        for (unsigned int i = 0; i < layers.size(); ++i)
        {
            delete layers[i];
        }
    }
    void simulate()
    {
        if (nudge_duration > 0.0)
        {
            nudge_duration -= SIM_TIME_PER_PHYSICS_TICK;
            if (nudge_duration <= 0.0)
            {
                nudge_duration = 0.0;
                nudge(-current_nudge.x, -current_nudge.y, 0.0);
            }
        }
        for (unsigned int i = 0; i < tilt_bobs.size(); ++i)
        {
            tilt_bobs[i]->simulate();
            pump_triggers(tilt_bobs[i]);
        }
        for (unsigned int i = 0; i < layers.size(); ++i)
        {
            layers[i]->simulate();
            pump_triggers(layers[i]);
        }
        for (unsigned int i = 0; i < spinners.size(); ++i)
        {
            spinners[i]->simulate();
            pump_triggers(spinners[i]);
        }
        if (backglass)
        {
            backglass->pass_time(MS_PER_FRAME / PHYSICS_TICKS_PER_FRAME);
            pump_triggers(backglass);
        }
        global_time += SIM_TIME_PER_PHYSICS_TICK;
    }
    void nudge(double vx, double vy, double duration)
    {
        Point nudge_vector = Point(vx, vy);
        Point nudge_origin = Point(dimensions.x / 2.0, dimensions.y);
        Point table_center_mass = Point(dimensions.x / 2.0, 0.0);
        if (nudge_duration > 0.0)
        {
            nudge_duration = 0.0;
            nudge(-current_nudge.x, -current_nudge.y, 0.0);
            return;
        }
        for (unsigned int i = 0; i < tilt_bobs.size(); ++i)
        {
            Point positional_nudge_vector = Point(0.0, 0.0);
            TiltBob* tilt_bob = tilt_bobs[i];
            convert_nudge(&nudge_vector, &(tilt_bob->position), &nudge_origin, &table_center_mass, &positional_nudge_vector);
            tilt_bob->velocity.x -= positional_nudge_vector.x;
            tilt_bob->velocity.y -= positional_nudge_vector.y;
        }
        for (unsigned int i = 0; i < spinners.size(); ++i)
        {
            //TODO (Phase X) test spinner nudging
            Point positional_nudge_vector = Point(0.0, 0.0);
            Spinner* spinner = spinners[i];
            convert_nudge(&nudge_vector, &(spinner->center_mass), &nudge_origin, &table_center_mass, &positional_nudge_vector);
            double aligned_nudge_magnitude = positional_nudge_vector.dot_product(spinner->unit_normal);
            double nudge_angular_velocity = aligned_nudge_magnitude / spinner->center_mass_radius * std::cos(spinner->angle);
            spinner->angular_velocity += nudge_angular_velocity;
        }
        for (unsigned int i = 0; i < layers.size(); ++i)
        {
            Layer* layer = layers[i];
            for (unsigned int j = 0; j < layer->balls.size(); ++j)
            {
                Point positional_nudge_vector = Point(0.0, 0.0);
                Ball* ball = layer->balls[j].first;
                if (!uniform_nudge)
                {
                    convert_nudge(&nudge_vector, &(ball->position), &nudge_origin, &table_center_mass, &positional_nudge_vector);
                } else
                {
                    positional_nudge_vector.x = nudge_vector.x;
                    positional_nudge_vector.y = nudge_vector.y;
                }
                //TODO rotate velocity vector of ball for horizontal nudges?
                ball->velocity.x -= positional_nudge_vector.x;
                ball->velocity.y -= positional_nudge_vector.y;
                //TODO (Phase X) extra testing on nudges affecting ball rotation
                Point table_center_mass_to_ball = Point(ball->position.x - table_center_mass.x, ball->position.y - table_center_mass.y);
                Point center_mass_to_nudge_origin = Point(nudge_origin.x - table_center_mass.x, nudge_origin.y - table_center_mass.y);
                ball->angular_velocity += nudge_vector.x * (table_center_mass_to_ball.magnitude() + ball->radius) / center_mass_to_nudge_origin.magnitude() / (ball->radius * 2.0 * PI);
            }
        }
        if (duration > 0.0)
        {
            current_nudge.x = vx;
            current_nudge.y = vy;
        }
        nudge_duration = duration;
    }
    void pump_triggers(Layer* layer)
    {
        while (layer->trigger_queue.size())
        {
            trigger_queue.push_back(layer->trigger_queue[0]);
            layer->trigger_queue.erase(layer->trigger_queue.begin());
        }
    }
    void pump_triggers(TiltBob* tilt_bob)
    {
        while (tilt_bob->trigger_queue.size())
        {
            trigger_queue.push_back(tilt_bob->trigger_queue[0]);
            tilt_bob->trigger_queue.erase(tilt_bob->trigger_queue.begin());
        }
    }
    void pump_triggers(Spinner* spinner)
    {
        while (spinner->trigger_queue.size())
        {
            trigger_queue.push_back(spinner->trigger_queue[0]);
            spinner->trigger_queue.erase(spinner->trigger_queue.begin());
        }
    }
    void pump_triggers(Backglass* backglass)
    {
        while (backglass->trigger_queue.size())
        {
            trigger_queue.push_back(backglass->trigger_queue.front());
            backglass->trigger_queue.pop_front();
        }
    }
    void render(SDL_Renderer* renderer, View* view)
    {
        for (unsigned int i = 0; i < layers.size(); ++i)
        {
            Layer* current_layer = layers[i];
            if (current_layer)
            {
                current_layer->render(renderer, view);
            }
        }
        if (view->debug)
        {
            for (unsigned int i = 0; i < tilt_bobs.size(); ++i)
            {
                TiltBob* current_bob = tilt_bobs[i];
                if (current_bob)
                {
                    current_bob->render(renderer, view);
                }
            }
            for (unsigned int i = 0; i < spinners.size(); ++i)
            {
                Spinner* current_spinner = spinners[i];
                if (current_spinner)
                {
                    current_spinner->render(renderer, view);
                }
            }
        }
        if (backglass)
        {
            backglass->render(renderer, view);
            backglass->HUD_render(renderer, view);
        }
    }
};

void trim_line(std::string& line)
{
    bool in_quotes = false;
    for (unsigned int i = 0; i < line.size(); ++i)
    {
        if (line[i] == '"' && (i == 0 || line[i - 1] != '\\'))
        {
            in_quotes = !in_quotes;
        }
        if (!in_quotes && line[i] == '\'')
        {
            line.resize(i);
        }
    }
    while (line.size() && (line.back() == '\n' || line.back() == ' ' || line.back() == '\t'))
    {
        line.pop_back();
    }
}

class Script;

void save_config(const char filename[], Script* script);
void load_config(const char filename[], Script* script);
void default_config(Script* script);

class Script
{
public:
    SDL_Renderer* renderer = nullptr;
    SDL_Window* window = nullptr;
    Syntax* syntax = nullptr;
    View* view = nullptr;
    SDL_Texture* background = nullptr;
    int bg_width = 0;
    int bg_height = 0;
    Table* table = nullptr;
    InputHandler* input_handler = nullptr;
    Mix_Music* song_playing = nullptr;
    int volume = 128;
    int music_volume = 128;
    int sound_volume = 128;
    bool mute = false;
    Point turtle_position = Point(0.0, 0.0);
    double turtle_heading = 0.0;
    unsigned int console_entries = 8;
    bool show_console = false;
    std::deque<std::string> console;
    std::string pending_console_line;
    std::vector<Instruction*> instructions;
    std::unordered_map<std::string, std::pair<DataType, int>> variable_index_table;
    std::vector<int> integer_variables;
    std::vector<double> double_variables;
    std::vector<std::string> string_variables;
    std::vector<std::pair<int, int*>> arrays; //pair is size, pointer
    std::vector<Collider*> colliders;
    std::vector<Sprite*> sprites;
    std::vector<Sound*> sounds;
    std::vector<Mix_Music*> songs;
    std::vector<Ball*> balls;
    std::vector<Magnet*> magnets;
    std::unordered_map<Event, std::unordered_map<int, int>*> trigger_points; //second map pairs parameter to script line
    std::vector<SpriteInstance*> sprite_instances; //Note: cleanup is handled by Layer
    std::vector<DotMatrixDisplay*> dot_matrix_displays;
    std::vector<AlphanumericDisplay*> alphanumeric_displays;
    std::vector<ScoreReel*> score_reels;
    std::vector<Backglass*> backglasses;
    std::vector<Font*> fonts;
    std::vector<FontSize*> font_sizes;
    std::vector<Record*> records;
    std::deque<std::pair<Event, int>> trigger_queue;
    std::unordered_map<int, bool> timers; //maps integer variable index to timer running bool
    int next_switch_id = 0;
    Script(SDL_Renderer* renderer, SDL_Window* window, Syntax* syntax, std::string filename) : renderer(renderer), syntax(syntax)
    {
        std::ifstream file(filename);
        std::string line;
        int line_number = 0;
        bool loaded = false;
        while (std::getline(file, line))
        {
            //TODO (Phase VIII) abort parsing on error?
            trim_line(line);
            if (line.size() == 0)
            {
                continue;
            }
            if (line[line.size() - 1] == ':')
            {
                parse_trigger(line, line_number);
            } else
            {
                if (parse_instruction(line) == 0)
                {
                    ++line_number;
                    if (!loaded && instructions.back()->command == BREAK)
                    {
                        add_trigger(LOAD, 0);
                        resolve_queue();
                        loaded = true;
                    }
                    if (instructions.back()->command == LABEL)
                    {
                        integer_variables[instructions.back()->int_args[0]] = line_number - 1;
                    }
                }
            }
        }
    }
    ~Script()
    {
        //TODO (Phase X) test against memory leaks
        if (background)
        {
            SDL_DestroyTexture(background);
        }
        for (unsigned int i = 0; i < instructions.size(); ++i)
        {
            delete instructions[i];
        }
        for (unsigned int i = 0; i < arrays.size(); ++i)
        {
            delete[] arrays[i].second;
        }
        for (unsigned int i = 0; i < colliders.size(); ++i)
        {
            delete colliders[i];
        }
        for (unsigned int i = 0; i < sprites.size(); ++i)
        {
            delete sprites[i];
        }
        for (unsigned int i = 0; i < sounds.size(); ++i)
        {
            delete sounds[i];
        }
        for (unsigned int i = 0; i < balls.size(); ++i)
        {
            delete balls[i];
        }
        for (unsigned int i = 0; i < magnets.size(); ++i)
        {
            delete magnets[i];
        }
        for (unsigned int i = 0; i < dot_matrix_displays.size(); ++i)
        {
            delete dot_matrix_displays[i];
        }
        for (unsigned int i = 0; i < alphanumeric_displays.size(); ++i)
        {
            delete alphanumeric_displays[i];
        }
        for (unsigned int i = 0; i < score_reels.size(); ++i)
        {
            delete score_reels[i];
        }
        for (unsigned int i = 0; i < backglasses.size(); ++i)
        {
            delete backglasses[i];
        }
        for (unsigned int i = 0; i < fonts.size(); ++i)
        {
            delete fonts[i];
        }
        for (unsigned int i = 0; i < font_sizes.size(); ++i)
        {
            delete font_sizes[i];
        }
        for (unsigned int i = 0; i < records.size(); ++i)
        {
            delete records[i];
        }
        for (auto i = trigger_points.begin(); i != trigger_points.end(); ++i)
        {
            delete i->second;
        }
        if (table)
        {
            delete table;
        }
    }
    void console_log(std::string text)
    {
        std::cout << text << std::endl;
        console.push_back(text);
        if (console.size() > console_entries)
        {
            console.pop_front();
        }
    }
    void process_event(SDL_Event* e)
    {
        switch (e->type)
        {
        case SDL_KEYDOWN:
            if (e->key.keysym.sym == SDLK_RETURN)
            {
                console.push_back(pending_console_line);
                if (pending_console_line.size())
                {
                    trim_line(pending_console_line);
                    parse_instruction(pending_console_line, true);
                }
                if (console.size() > console_entries)
                {
                    console.pop_front();
                }
                pending_console_line.clear();
            } else if (e->key.keysym.sym == SDLK_BACKSPACE)
            {
                if (pending_console_line.size())
                {
                    pending_console_line.pop_back();
                }
            } else if(e->key.keysym.sym == SDLK_x && SDL_GetModState() & KMOD_CTRL)
            {
                SDL_SetClipboardText(pending_console_line.c_str());
                pending_console_line.clear();
            } else if(e->key.keysym.sym == SDLK_c && SDL_GetModState() & KMOD_CTRL)
            {
                SDL_SetClipboardText(pending_console_line.c_str());
            } else if(e->key.keysym.sym == SDLK_v && SDL_GetModState() & KMOD_CTRL)
            {
                char* temp = SDL_GetClipboardText();
                std::string temp_string = temp;
                for (unsigned int i = 0; i < temp_string.size(); ++i)
                {
                    if (temp_string[i] < 32 || temp_string[i] > 126)
                    {
                        temp_string.erase(temp_string.begin() + i);
                        --i;
                    }
                }
                pending_console_line.append(temp_string.c_str());
                SDL_free(temp);
            }
            break;
        case SDL_TEXTINPUT:
            if (e && strcmp(e->text.text, "`"))
            {
                std::string temp_string = e->text.text;
                for (unsigned int i = 0; i < temp_string.size(); ++i)
                {
                    if (temp_string[i] < 32 || temp_string[i] > 126)
                    {
                        temp_string.erase(temp_string.begin() + i);
                        --i;
                    }
                }
                pending_console_line.append(temp_string);
            }
            break;
        case SDL_TEXTEDITING:
            break;
        }
    }
    void render_console(SDL_Renderer* renderer, Font* font, int x, int y, int line_height)
    {
        //Note: line height is applied before UI scale
        for (unsigned int i = 0; i < console.size(); ++i)
        {
            font->show_text(renderer, console[i], x, y + line_height * i, ALIGN_TOP_LEFT, DISPLAY_WIDTH, line_height, UI_SCALE);
        }
        std::string display_string = ">";
        display_string.append(pending_console_line);
        font->show_text(renderer, display_string, x, y + line_height * UI_SCALE * console.size(), ALIGN_TOP_LEFT, DISPLAY_WIDTH, line_height, UI_SCALE);
    }
    void refresh_instruction(Instruction* instruction)
    {
        for (unsigned int i = 0; i < instruction->variable_table.size(); ++i)
        {
            auto variable_info = instruction->variable_table[i];
            switch (std::get<0>(variable_info))
            {
            case TYPE_INTEGER_VARIABLE:
                instruction->int_args[std::get<2>(variable_info)] = integer_variables[std::get<1>(variable_info)];
                break;
            case TYPE_DOUBLE_VARIABLE:
                instruction->double_args[std::get<2>(variable_info)] = double_variables[std::get<1>(variable_info)];
                break;
            case TYPE_STRING_VARIABLE:
                instruction->string_args[std::get<2>(variable_info)] = string_variables[std::get<1>(variable_info)];
                break;
            default:
                break;
            }
        }
    }
    void add_trigger(Event event, int parameter)
    {
        trigger_queue.push_back(std::pair<Event, int>(event, parameter));
    }
    void parse_trigger(std::string line, int line_number)
    {
        if (line.size() == 0)
        {
            return;
        }
        std::string word;
        std::istringstream linestream(line);
        linestream >> word;
        if (word.back() == ':')
        {
            word.pop_back();
            if (syntax->event_table.find(word) == syntax->event_table.end())
            {
                console_log("SCRIPT PARSING ERROR: undefined event");
                std::string error_message = "Line: ";
                error_message.append(line);
                console_log(error_message);
            } else
            {
                Event event = syntax->event_table[word];
                if (event == LOAD && line_number != 0)
                {
                    console_log("SCRIPT PARSING ERROR: LOAD trigger not on line 0");
                } else if (trigger_points.find(event) == trigger_points.end())
                {
                    std::unordered_map<int, int>* secondary_map = new std::unordered_map<int, int>;
                    (*secondary_map)[0] = line_number;
                    trigger_points[event] = secondary_map;
                } else
                {
                    console_log("SCRIPT PARSING ERROR: duplicate trigger");
                    std::string error_message = "Line: ";
                    error_message.append(line);
                    console_log(error_message);
                }
            }
        } else
        {
            Event event = syntax->event_table[word];
            linestream >> word;
            word.pop_back();
            int parameter;
            if (syntax->enum_table.find(word) != syntax->enum_table.end())
            {
                //TODO handle error when enum is not found in place of integer (currently fails in stoi and aborts program)
                parameter = syntax->enum_table[word];
            } else if (variable_index_table.find(word) != variable_index_table.end() && variable_index_table[word].first == TYPE_INTEGER_VARIABLE)
            {
                if (event == TIMER_DONE)
                {
                    if (variable_index_table[word].first == TYPE_INTEGER_VARIABLE && timers.find(variable_index_table[word].second) != timers.end())
                    {
                        //TODO (Phase X) test TIMER_DONE trigger
                        parameter = variable_index_table[word].second;
                    } else
                    {
                        console_log("SCRIPT PARSING ERROR: non-timer variable passed to TIMER_DONE trigger");
                        std::string error_message = "Line: ";
                        error_message.append(line);
                        console_log(error_message);
                    }
                } else
                {
                    parameter = integer_variables[variable_index_table[word].second];
                }
            } else
            {
                parameter = stoi(word);
            }
            if (trigger_points.find(event) == trigger_points.end())
            {
                std::unordered_map<int, int>* secondary_map = new std::unordered_map<int, int>;
                (*secondary_map)[parameter] = line_number;
                trigger_points[event] = secondary_map;
            } else
            {
                std::unordered_map<int, int>* secondary_map = trigger_points[event];
                if (secondary_map->find(parameter) == secondary_map->end())
                {
                    (*secondary_map)[parameter] = line_number;
                } else
                {
                    console_log("SCRIPT PARSING ERROR: duplicate trigger");
                    std::string error_message = "Line: ";
                    error_message.append(line);
                    console_log(error_message);
                }
            }
        }
    }
    int parse_instruction(std::string line, bool instant_run = false)
    {
        if (line.size() == 0)
        {
            //Line contains only comments/whitespace
            return 0;
        }
        Instruction* instruction = new Instruction;
        std::string word;
        std::vector<std::string> words;
        std::istringstream linestream(line);
        while (linestream >> word)
        {
            while (word[0] == '"' && (!(word[word.size() - 1] == '"') || (word.size() > 1 && word[word.size() - 2] == '\\')))
            {
                std::string next;
                if (linestream >> next)
                {
                    word.append(" ");
                    word.append(next);
                } else
                {
                    console_log("SCRIPT PARSING ERROR: string contains end of line");
                    std::string error_message = "Line: ";
                    error_message.append(line);
                    console_log(error_message);
                    delete instruction;
                    return -1;
                }
            }
            words.push_back(word);
        }
        if (syntax->command_table.find(words[0]) == syntax->command_table.end())
        {
            std::string error_message = "SCRIPT PARSING ERROR: ";
            error_message.append(words[0]);
            error_message.append(" is not a command");
            console_log(error_message);
            error_message = "Line: ";
            error_message.append(line);
            console_log(error_message);
            delete instruction;
            return -1;
        } else
        {
            instruction->command = syntax->command_table[words[0]];
            std::vector<DataType>* args = syntax->arg_table[instruction->command];
            if (words.size() != args->size() + 1)
            {
                std::string error_message = "SCRIPT PARSING ERROR: ";
                error_message.append(std::to_string(words.size() - 1));
                error_message.append(" argument(s) provided; ");
                error_message.append(std::to_string(args->size()));
                error_message.append(" expected");
                console_log(error_message);
                error_message = "Line: ";
                error_message.append(line);
                console_log(error_message);
                delete instruction;
                return -1;
            }
            for (unsigned int i = 0; i < args->size(); ++i)
            {
                std::string arg_string = words[i + 1];
                DataType arg_type;
                if (arg_string[0] == '"')
                {
                    arg_type = TYPE_STRING;
                } else if ((arg_string[0] >= '0' && arg_string[0] <= '9') || arg_string[0] == '.' || arg_string[0] == '-' || syntax->enum_table.find(arg_string) != syntax->enum_table.end())
                {
                    if (arg_string.find('.') != std::string::npos)
                    {
                        arg_type = TYPE_DOUBLE;
                    } else
                    {
                        arg_type = TYPE_INTEGER;
                    }
                } else
                {
                    auto arg_info = variable_index_table.find(arg_string);
                    if (arg_info != variable_index_table.end())
                    {
                        arg_type = arg_info->second.first;
                    } else
                    {
                        std::string error_message = "SCRIPT PARSING ERROR: ";
                        error_message.append(arg_string);
                        error_message.append(" not defined");
                        console_log(error_message);
                        error_message = "Line: ";
                        error_message.append(line);
                        console_log(error_message);
                        delete instruction;
                        return -1;
                    }
                }
                switch ((*args)[i])
                {
                case TYPE_INTEGER:
                    if (arg_type == TYPE_INTEGER)
                    {
                        if (syntax->enum_table.find(arg_string) != syntax->enum_table.end())
                        {
                            instruction->int_args.push_back(syntax->enum_table[arg_string]);
                        } else
                        {
                            try
                            {
                                instruction->int_args.push_back(std::stoi(arg_string));
                            } catch (std::invalid_argument &error)
                            {
                                console_log("SCRIPT PARSING ERROR: argument cannot be parsed");
                                std::string error_message = "Line: ";
                                error_message.append(line);
                                console_log(error_message);
                                delete instruction;
                                return -1;
                            } catch (std::out_of_range &error)
                            {
                                console_log("SCRIPT PARSING ERROR: argument out of range");
                                std::string error_message = "Line: ";
                                error_message.append(line);
                                console_log(error_message);
                                delete instruction;
                                return -1;
                            }
                        }
                    } else if (arg_type == TYPE_INTEGER_VARIABLE)
                    {
                        instruction->int_args.push_back(integer_variables[variable_index_table.find(arg_string)->second.second]);
                        instruction->variable_table.push_back(std::tuple<DataType, int, int>(TYPE_INTEGER_VARIABLE, variable_index_table[arg_string].second, instruction->int_args.size() - 1));
                    } else
                    {
                        console_log("SCRIPT PARSING ERROR: type mismatch (expecting integer)");
                        std::string error_message = "Line: ";
                        error_message.append(line);
                        console_log(error_message);
                        delete instruction;
                        return -1;
                    }
                    break;
                case TYPE_DOUBLE:
                    if (arg_type == TYPE_DOUBLE)
                    {
                        //TODO allow angles in degrees (converted to radians on parse)
                        try
                        {
                            instruction->double_args.push_back(std::stod(arg_string));
                        } catch (std::invalid_argument &error)
                        {
                            console_log("SCRIPT PARSING ERROR: argument cannot be parsed");
                            std::string error_message = "Line: ";
                            error_message.append(line);
                            console_log(error_message);
                            delete instruction;
                            return -1;
                        } catch (std::out_of_range &error)
                        {
                            console_log("SCRIPT PARSING ERROR: argument out of range");
                            std::string error_message = "Line: ";
                            error_message.append(line);
                            console_log(error_message);
                            delete instruction;
                            return -1;
                        }
                    } else if (arg_type == TYPE_DOUBLE_VARIABLE)
                    {
                        instruction->double_args.push_back(double_variables[variable_index_table.find(arg_string)->second.second]);
                        instruction->variable_table.push_back(std::tuple<DataType, int, int>(TYPE_DOUBLE_VARIABLE, variable_index_table[arg_string].second, instruction->double_args.size() - 1));
                    } else
                    {
                        console_log("SCRIPT PARSING ERROR: type mismatch (expecting double)");
                        std::string error_message = "Line: ";
                        error_message.append(line);
                        console_log(error_message);
                        delete instruction;
                        return -1;
                    }
                    break;
                case TYPE_STRING:
                    if (arg_type == TYPE_STRING)
                    {
                        if (arg_string[0] == '"')
                        {
                            arg_string.erase(arg_string.begin());
                        }
                        for (unsigned int i = 0; i < arg_string.size() - 1; ++i)
                        {
                            if (arg_string[i] == '\\' && arg_string[i + 1] == '"')
                            {
                                arg_string.erase(arg_string.begin() + i);
                            }
                        }
                        if (arg_string[arg_string.size() - 1] == '"')
                        {
                            arg_string.pop_back();
                        }
                        instruction->string_args.push_back(arg_string);
                    } else if (arg_type == TYPE_STRING_VARIABLE)
                    {
                        instruction->string_args.push_back(string_variables[variable_index_table.find(arg_string)->second.second]);
                        instruction->variable_table.push_back(std::tuple<DataType, int, int>(TYPE_STRING_VARIABLE, variable_index_table[arg_string].second, instruction->string_args.size() - 1));
                    } else
                    {
                        console_log("SCRIPT PARSING ERROR: type mismatch (expecting string)");
                        std::string error_message = "Line: ";
                        error_message.append(line);
                        console_log(error_message);
                        delete instruction;
                        return -1;
                    }
                    break;
                case TYPE_INTEGER_VARIABLE:
                    if (arg_type == TYPE_INTEGER_VARIABLE)
                    {
                        instruction->int_args.push_back(variable_index_table[arg_string].second);
                    } else if (arg_type == TYPE_INTEGER && (arg_string == "-1" || arg_string == "NULL"))
                    {
                        //TODO refactor or remove NULL for all variable types (currently doesn't parse)
                        instruction->int_args.push_back(-1);
                    } else
                    {
                        console_log("SCRIPT PARSING ERROR: type mismatch (expecting integer variable)");
                        std::string error_message = "Line: ";
                        error_message.append(line);
                        console_log(error_message);
                        delete instruction;
                        return -1;
                    }
                    break;
                case TYPE_DOUBLE_VARIABLE:
                    if (arg_type == TYPE_DOUBLE_VARIABLE)
                    {
                        instruction->int_args.push_back(variable_index_table[arg_string].second);
                    } else if (arg_type == TYPE_INTEGER && (arg_string == "-1" || arg_string == "NULL"))
                    {
                        instruction->int_args.push_back(-1);
                    } else
                    {
                        console_log("SCRIPT PARSING ERROR: type mismatch (expecting double variable)");
                        std::string error_message = "Line: ";
                        error_message.append(line);
                        console_log(error_message);
                        delete instruction;
                        return -1;
                    }
                    break;
                case TYPE_STRING_VARIABLE:
                    if (arg_type == TYPE_STRING_VARIABLE)
                    {
                        instruction->int_args.push_back(variable_index_table[arg_string].second);
                    } else if (arg_type == TYPE_INTEGER && (arg_string == "-1" || arg_string == "NULL"))
                    {
                        instruction->int_args.push_back(-1);
                    } else
                    {
                        console_log("SCRIPT PARSING ERROR: type mismatch (expecting string variable)");
                        std::string error_message = "Line: ";
                        error_message.append(line);
                        console_log(error_message);
                        delete instruction;
                        return -1;
                    }
                    break;
                case TYPE_ARRAY:
                    if (arg_type == TYPE_ARRAY)
                    {
                        instruction->int_args.push_back(variable_index_table[arg_string].second);
                    } else
                    {
                        console_log("SCRIPT PARSING ERROR: type mismatch (expecting array)");
                        std::string error_message = "Line: ";
                        error_message.append(line);
                        console_log(error_message);
                        delete instruction;
                        return -1;
                    }
                    break;
                default:
                    break;
                }
            }
            if (instant_run || instruction->command == INTEGER || instruction->command == DOUBLE || instruction->command == STRING || instruction->command == ARRAY)
            {
                execute(instruction);
            }
            if (instant_run)
            {
                delete instruction;
            } else
            {
                instructions.push_back(instruction);
            }
            return 0;
        }
    }
    void pass_balls()
    {
        for (unsigned int i = 0; i < balls.size(); ++i)
        {
            Ball* ball = balls[i];
            if (ball->layer != ball->layer_pending)
            {
                table->layers[ball->layer]->pass_ball(ball, table->layers[ball->layer_pending]);
            }
        }
    }
    int execute(Instruction* instruction)
    {
        //TODO (Phase VIII) wherever possible, provide default arguments when -1 is passed
        //[IMPLEMENTATION LIST]
        refresh_instruction(instruction);
        switch (instruction->command)
        {
        //TODO (Phase X) double-check all commands for bounds checking, proper argument passing, and no memory leaks
        case DEBUG_RENDERER:
            if (view)
            {
                view->debug = instruction->int_args[0];
            }
            break;
        case SPRITE_RENDERER:
            if (view)
            {
                view->sprites = instruction->int_args[0];
            }
            break;
        case HUD_RENDERER:
            if (view)
            {
                view->HUD = instruction->int_args[0];
            }
            break;
        case SET_VIEW:
            if (view)
            {
                view->pending_center.x = instruction->double_args[0];
                view->pending_center.y = instruction->double_args[1];
                view->pending_zoom = instruction->double_args[2];
                view->transition_time = instruction->double_args[3];
            }
            break;
        case GET_VIEW:
            if (view)
            {
                int x_index = instruction->int_args[0];
                if (x_index != -1)
                {
                    double_variables[x_index] = view->center.x;
                }
                int y_index = instruction->int_args[1];
                if (y_index != -1)
                {
                    double_variables[y_index] = view->center.y;
                }
                int zoom_index = instruction->int_args[2];
                if (zoom_index != -1)
                {
                    double_variables[zoom_index] = view->zoom;
                }
            }
            break;
        case TABLE:
            if (!table)
            {
                table = new Table;
                table->dimensions.x = instruction->double_args[0];
                table->dimensions.y = instruction->double_args[1];
                table->uniform_nudge = instruction->int_args[0];
            } else
            {
                console_log("SCRIPT ERROR: table already declared");
                return -1;
            }
            break;
        case LAYER:
            {
                Layer* layer = new Layer;
                layer->top_left.x = instruction->double_args[0];
                layer->top_left.y = instruction->double_args[1];
                layer->dimensions.x = instruction->double_args[2];
                layer->dimensions.y = instruction->double_args[3];
                layer->gradient_mode = static_cast<GradientMode>(instruction->int_args[0]);
                layer->gradient_vector.x = instruction->double_args[4];
                layer->gradient_vector.y = instruction->double_args[5];
                layer->center.x = instruction->double_args[6];
                layer->center.y = instruction->double_args[7];
                if (table)
                {
                    table->layers.push_back(layer);
                    layer->id = table->layers.size() - 1;
                    int index = instruction->int_args[1];
                    if (index != -1)
                    {
                        integer_variables[index] = table->layers.size() - 1;
                    }
                } else
                {
                    delete layer;
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
            }
            break;
        case POINT:
            {
                colliders.push_back(new PointCollider(instruction->double_args[0], instruction->double_args[1]));
                colliders.back()->id = colliders.size() - 1;
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = colliders.size() - 1;
                }
            }
            break;
        case LINE:
            {
                colliders.push_back(new LineCollider(instruction->double_args[0], instruction->double_args[1], instruction->double_args[2], instruction->double_args[3]));
                colliders.back()->id = colliders.size() - 1;
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = colliders.size() - 1;
                }
            }
            break;
        case ARC:
            {
                colliders.push_back(new ArcCollider(instruction->double_args[0], instruction->double_args[1], instruction->double_args[2], instruction->double_args[3], instruction->double_args[4]));
                colliders.back()->id = colliders.size() - 1;
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = colliders.size() - 1;
                }
            }
            break;
        case CIRCLE:
            {
                colliders.push_back(new ArcCollider(instruction->double_args[0], instruction->double_args[1], instruction->double_args[2]));
                colliders.back()->id = colliders.size() - 1;
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = colliders.size() - 1;
                }
            }
            break;
        case TURTLE_TELEPORT:
            {
                turtle_position.x = instruction->double_args[0];
                turtle_position.y = instruction->double_args[1];
            }
            break;
        case TURTLE_SET_HEADING:
            {
                turtle_heading = standardize_angle(instruction->double_args[0]);
            }
            break;
        case TURTLE_FORWARD:
            {
                double distance = instruction->double_args[0];
                turtle_position.x += distance * std::cos(turtle_heading);
                turtle_position.y += distance * std::sin(turtle_heading);
            }
            break;
        case TURTLE_TURN:
            {
                turtle_heading += instruction->double_args[0];
            }
            break;
        case TURTLE_POINT:
            {
                colliders.push_back(new PointCollider(turtle_position.x, turtle_position.y));
                colliders.back()->id = colliders.size() - 1;
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = colliders.size() - 1;
                }
            }
            break;
        case TURTLE_LINE:
            {
                Point line_start = Point(turtle_position.x, turtle_position.y);
                double distance = instruction->double_args[0];
                turtle_position.x += distance * std::cos(turtle_heading);
                turtle_position.y += distance * std::sin(turtle_heading);
                colliders.push_back(new LineCollider(line_start.x, line_start.y, turtle_position.x, turtle_position.y));
                colliders.back()->id = colliders.size() - 1;
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = colliders.size() - 1;
                }
            }
            break;
        case TURTLE_ARC:
            {
                //TODO (Phase X) test when turtle heading or arc width is outside [-2pi, 2pi]
                double radius = instruction->double_args[0];
                double angle = instruction->double_args[1];
                Point center = Point(turtle_position.x, turtle_position.y);
                double arc_start_angle = standardize_angle(turtle_heading);
                if (angle > 0.0)
                {
                    center.x += radius * std::cos(turtle_heading + PI / 2.0);
                    center.y += radius * std::sin(turtle_heading + PI / 2.0);
                    arc_start_angle -= PI / 2.0;
                } else
                {
                    center.x += radius * std::cos(turtle_heading - PI / 2.0);
                    center.y += radius * std::sin(turtle_heading - PI / 2.0);
                    arc_start_angle += PI / 2.0;
                }
                double arc_end_angle = arc_start_angle + angle;
                if (arc_end_angle > 2 * PI)
                {
                    arc_start_angle -= 2 * PI;
                    arc_end_angle -= 2 * PI;
                }
                if (arc_end_angle > arc_start_angle)
                {
                    colliders.push_back(new ArcCollider(center.x, center.y, radius, arc_start_angle, arc_end_angle));
                } else
                {
                    colliders.push_back(new ArcCollider(center.x, center.y, radius, arc_end_angle, arc_start_angle));
                }
                colliders.back()->id = colliders.size() - 1;
                turtle_position.x = center.x + radius * std::cos(arc_end_angle);
                turtle_position.y = center.y + radius * std::sin(arc_end_angle);
                turtle_heading = standardize_angle(turtle_heading + angle);
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = colliders.size() - 1;
                }
            }
            break;
        case TURTLE_GET_POSITION:
            {
                int x_index = instruction->int_args[0];
                if (x_index != -1)
                {
                    double_variables[x_index] = turtle_position.x;
                }
                int y_index = instruction->int_args[1];
                if (y_index != -1)
                {
                    double_variables[y_index] = turtle_position.y;
                }
            }
            break;
        case BALL:
            {
                if (table)
                {
                    unsigned int layer = instruction->int_args[0];
                    if (layer >= table->layers.size())
                    {
                        console_log("SCRIPT ERROR: layer out of range");
                        return -1;
                    }
                    Ball* ball = new Ball(instruction->double_args[0], instruction->double_args[1]);
                    double radius = instruction->double_args[2];
                    if (radius > 0)
                    {
                        ball->radius = radius;
                    }
                    double mass = instruction->double_args[3];
                    if (mass > 0)
                    {
                        ball->mass = mass;
                    }
                    double moment_ratio = instruction->double_args[4];
                    if (moment_ratio <= 0.0 || moment_ratio > 1.0)
                    {
                        moment_ratio = 0.4;
                    }
                    ball->moment_ratio = moment_ratio;
                    int magnetic = instruction->int_args[1];
                    ball->magnetic = magnetic;
                    bool show_spin = instruction->int_args[2];
                    ball->show_spin = show_spin;
                    int sprite_index = instruction->int_args[3];
                    Sprite* sprite = nullptr;
                    if (sprite_index >= 0 && static_cast<unsigned int>(sprite_index) < sprites.size())
                    {
                        sprite = sprites[sprite_index];
                    } else
                    {
                        if (sprite_index != -1)
                        {
                            console_log("SCRIPT ERROR: spritesheet out of range");
                            return -1;
                        }
                    }
                    int impact_sound_index = instruction->int_args[4];
                    if (impact_sound_index >= 0 && static_cast<unsigned int>(impact_sound_index) < sounds.size())
                    {
                        ball->impact_sound = sounds[impact_sound_index];
                    } else
                    {
                        if (impact_sound_index != -1)
                        {
                            console_log("SCRIPT ERROR: sound out of range");
                            return -1;
                        }
                    }
                    int roll_sound_index = instruction->int_args[5];
                    if (roll_sound_index >= 0 && static_cast<unsigned int>(roll_sound_index) < sounds.size())
                    {
                        ball->roll_sound = sounds[roll_sound_index];
                    } else
                    {
                        if (roll_sound_index != -1)
                        {
                            console_log("SCRIPT ERROR: sound out of range");
                            return -1;
                        }
                    }
                    ball->layer = layer;
                    table->layers[layer]->add_ball(ball, sprite);
                    balls.push_back(ball);
                    int index = instruction->int_args[6];
                    if (index != -1)
                    {
                        integer_variables[index] = balls.size() - 1;
                    }
                } else
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
            }
            break;
        case STANDARD_BALL: //layer, x, y, ID out
            if (table)
            {
                unsigned int layer = instruction->int_args[0];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                Ball* ball = new Ball(instruction->double_args[0], instruction->double_args[1]);
                //TODO (Phase X) test validity of reusing sprite for multiple standard balls; do so if possible
                //TODO (Phase VII) standard sprite
                //TODO (Phase VII) standard impact sound
                //TODO (Phase VII) standard roll sound
                ball->layer = layer;
                table->layers[layer]->add_ball(ball, nullptr);
                balls.push_back(ball);
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    integer_variables[index] = balls.size() - 1;
                }
            } else
            {
                console_log("SCRIPT ERROR: no table declared");
                return -1;
            }
            break;
        case MAGNET:
            {
                if (table)
                {
                    unsigned int layer = instruction->int_args[0];
                    if (layer >= table->layers.size())
                    {
                        console_log("SCRIPT ERROR: layer out of range");
                        return -1;
                    }
                    Magnet* magnet = new Magnet(0.0, 0.0, 0.0);
                    magnet->position.x = instruction->double_args[0];
                    magnet->position.y = instruction->double_args[1];
                    magnet->magnetism = instruction->double_args[2];
                    magnet->radius = instruction->double_args[3];
                    magnets.push_back(magnet);
                    table->layers[layer]->add_magnet(magnet);
                    double depth = instruction->double_args[4];
                    if (depth > 0.0)
                    {
                        magnet->depth = depth;
                    }
                    magnet->hole = instruction->int_args[1];
                    int index = instruction->int_args[2];
                    if (index != -1)
                    {
                        integer_variables[index] = magnets.size() - 1;
                    }
                } else
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
            }
            break;
        case SPRITE:
            {
                if (table)
                {
                    unsigned int layer = instruction->int_args[1];
                    if (layer >= table->layers.size())
                    {
                        console_log("SCRIPT ERROR: layer out of range");
                        return -1;
                    }
                    Sprite* sprite = nullptr;
                    if (instruction->int_args[0] >= 0 && static_cast<unsigned int>(instruction->int_args[0]) < sprites.size())
                    {
                        sprite = sprites[instruction->int_args[0]];
                    } else
                    {
                        console_log("SCRIPT ERROR: spritesheet out of range");
                        return -1;
                    }
                    SpriteInstance* sprite_instance = table->layers[layer]->add_sprite(sprite, instruction->double_args[0], instruction->double_args[1], instruction->int_args[2], 0.0);
                    sprite_instances.push_back(sprite_instance);
                    int index = instruction->int_args[3];
                    if (index != -1)
                    {
                        integer_variables[index] = sprite_instances.size() - 1;
                    }
                } else
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
            }
            break;
        case DMD_SPRITE: //ID in, DMD ID, x, y, ID out
            {
                int sprite_id = instruction->int_args[0];
                int dmd_id = instruction->int_args[1];
                int x = instruction->int_args[2];
                int y = instruction->int_args[3];
                int index = instruction->int_args[4];
                Sprite* sprite = nullptr;
                DotMatrixDisplay* dmd = nullptr;
                if (sprite_id >= 0 && static_cast<unsigned int>(sprite_id) < sprites.size())
                {
                    sprite = sprites[sprite_id];
                } else
                {
                    console_log("SCRIPT ERROR: spritesheet out of range");
                    return -1;
                }
                if (dmd_id >= 0 && static_cast<unsigned int>(dmd_id) < dot_matrix_displays.size())
                {
                    dmd = dot_matrix_displays[dmd_id];
                } else
                {
                    console_log("SCRIPT ERROR: dot matrix display out of range");
                    return -1;
                }
                SpriteInstance* sprite_instance = dmd->add_sprite(sprite, x, y);
                sprite_instances.push_back(sprite_instance);
                if (index != -1)
                {
                    integer_variables[index] = sprite_instances.size() - 1;
                }
            }
            break;
        case BACKGROUND:
            {
                if (background)
                {
                    SDL_DestroyTexture(background);
                    background = nullptr;
                }
                std::string filename = instruction->string_args[0];
                if (filename != "")
                {
                    std::ifstream file_test(filename.c_str());
                    if (!file_test.good())
                    {
                        std::string error_message = "SCRIPT ERROR: can't read from ";
                        error_message.append(filename);
                        console_log(error_message);
                        return -1;
                    }
                    SDL_Surface* background_surface = IMG_Load(filename.c_str());
                    background = SDL_CreateTextureFromSurface(renderer, background_surface);
                    SDL_FreeSurface(background_surface);
                    SDL_QueryTexture(background, nullptr, nullptr, &(bg_width), &(bg_height));
                }
            }
            break;
        case SPINNER: //collider ID in, layer, contact point radius, center mass radius, damping coefficient, gravity, sprite ID in, frames per rotation, ID out
            {
                if (table)
                {
                    Spinner* spinner = new Spinner;
                    if (instruction->int_args[0] >= 0 && static_cast<unsigned int>(instruction->int_args[0]) < colliders.size())
                    {
                        Collider* collider = colliders[instruction->int_args[0]];
                        LineCollider* line_collider = dynamic_cast<LineCollider*>(collider);
                        if (line_collider == nullptr)
                        {
                            console_log("SCRIPT ERROR: tried to make spinner from non-line collider");
                            delete spinner;
                            return -1;
                        } else
                        {
                            line_collider->spinner = spinner;
                            line_collider->solid = false;
                            spinner->center_mass.x = (line_collider->a.x + line_collider->b.x) / 2.0;
                            spinner->center_mass.y = (line_collider->a.y + line_collider->b.y) / 2.0;
                            Point a_to_b = Point(line_collider->b.x - line_collider->a.x, line_collider->b.y - line_collider->a.y);
                            Point a_to_b_unit = Point(a_to_b.x / a_to_b.magnitude(), a_to_b.y / a_to_b.magnitude());
                            spinner->unit_normal.x = -a_to_b_unit.y;
                            spinner->unit_normal.y = a_to_b_unit.x;
                        }
                        int layer = instruction->int_args[1];
                        if (layer >= 0 && static_cast<unsigned int>(layer) < table->layers.size())
                        {
                            table->layers[layer]->add_collider(line_collider);
                        } else
                        {
                            console_log("SCRIPT ERROR: layer out of range");
                            delete spinner;
                            return -1;
                        }
                    }
                    if (instruction->double_args[0] > 0.0)
                    {
                        spinner->contact_point_radius = instruction->double_args[0];
                    }
                    if (instruction->double_args[1] > 0.0)
                    {
                        spinner->center_mass_radius = instruction->double_args[1];
                    }
                    if (instruction->double_args[2] >= 0.0)
                    {
                        spinner->damping_coefficient = instruction->double_args[2];
                    }
                    if (instruction->double_args[3] >= 0.0)
                    {
                        spinner->gravity = instruction->double_args[3];
                    }
                    if (instruction->int_args[2] >= 0 && static_cast<unsigned int>(instruction->int_args[2]) < sprite_instances.size())
                    {
                        spinner->synced_sprite = sprite_instances[instruction->int_args[2]];
                        if (instruction->int_args[3] > 0)
                        {
                            spinner->frames_per_rotation = instruction->int_args[3];
                        }
                    } else if (instruction->int_args[2] >= 0)
                    {
                        console_log("SCRIPT ERROR: sprite out of range");
                        return -1;
                    }
                    spinner->id = table->spinners.size();
                    int index = instruction->int_args[4];
                    if (index != -1)
                    {
                        integer_variables[instruction->int_args[4]] = table->spinners.size();
                    }
                    table->spinners.push_back(spinner);
                } else
                {
                    console_log("SCRIPT ERROR: no table declared");
                }
            }
            break;
        case TILT_BOB: //x, y, tolerance, damping coefficient, bounce coefficient, gravity
            {
                if (table)
                {
                    TiltBob* tilt_bob = new TiltBob;
                    tilt_bob->position.x = instruction->double_args[0];
                    tilt_bob->position.y = instruction->double_args[1];
                    if (instruction->double_args[2] > 0.0)
                    {
                        tilt_bob->tolerance = instruction->double_args[2];
                    }
                    if (instruction->double_args[3] >= 0.0)
                    {
                        tilt_bob->damping_coefficient = instruction->double_args[3];
                    }
                    if (instruction->double_args[4] >= 1.0)
                    {
                        tilt_bob->bounce_coefficient = instruction->double_args[4];
                    }
                    if (instruction->double_args[5] >= 0.0)
                    {
                        tilt_bob->gravity = instruction->double_args[5];
                    }
                    tilt_bob->id = table->tilt_bobs.size();
                    table->tilt_bobs.push_back(tilt_bob);
                } else
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
            }
            break;
        case RESET_TILT_BOB:
            {
                if (table)
                {
                    TiltBob* tilt_bob = table->tilt_bobs[instruction->int_args[0]];
                    tilt_bob->bob_position.x = 0.0;
                    tilt_bob->bob_position.y = 0.0;
                    tilt_bob->velocity.x = 0.0;
                    tilt_bob->velocity.y = 0.0;
                    tilt_bob->trigger_queue.clear();
                } else
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
            }
            break;
        case SET_DETECTION_NORMAL:
            {
                int collider_index = instruction->int_args[0];
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    if (instruction->int_args[0] == -1)
                    {
                        collider_index = colliders.size() - 1;
                    } else
                    {
                        console_log("SCRIPT ERROR: collider out of range");
                        return -1;
                    }
                }
                Collider* collider = colliders[collider_index];
                double angle = instruction->double_args[0];
                Point detection_unit_normal = Point(std::cos(angle), std::sin(angle));
                collider->detection_unit_normal = detection_unit_normal;
                collider->has_detection_unit_normal = true;
            }
            break;
        case GRADIENT_VECTOR:
            if (table)
            {
                unsigned int layer = instruction->int_args[0];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                double x = instruction->double_args[0];
                double y = instruction->double_args[1];
                table->layers[layer]->gradient_vector.x = x;
                table->layers[layer]->gradient_vector.y = y;
            } else
            {
                console_log("SCRIPT ERROR: no table declared");
                return -1;
            }
            break;
        case BALL_RETURN:
            if (table)
            {
                unsigned int layer = instruction->int_args[0];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                double x = instruction->double_args[0];
                double y = instruction->double_args[1];
                table->layers[layer]->ball_return.x = x;
                table->layers[layer]->ball_return.y = y;
            } else
            {
                console_log("SCRIPT ERROR: no table declared");
                return -1;
            }
            break;
        case COMPOSITE:
            {
                if (instruction->int_args[1] < instruction->int_args[0])
                {
                    console_log("SCRIPT ERROR: composite collider range empty");
                    return -1;
                }
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[1]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                CompositeCollider* composite_collider = new CompositeCollider;
                for (int i = instruction->int_args[0]; i <= instruction->int_args[1]; ++i)
                {
                    Collider* collider = colliders[i];
                    if (dynamic_cast<PointCollider*>(collider) != nullptr)
                    {
                        composite_collider->add_collider(dynamic_cast<PointCollider*>(collider));
                    } else if (dynamic_cast<LineCollider*>(collider) != nullptr)
                    {
                        composite_collider->add_collider(dynamic_cast<LineCollider*>(collider));
                    } else if (dynamic_cast<ArcCollider*>(collider) != nullptr)
                    {
                        composite_collider->add_collider(dynamic_cast<ArcCollider*>(collider));
                    } else
                    {
                        console_log("SCRIPT ERROR: tried to add non-elementary collider to composite collider");
                        return -1;
                    }
                }
                composite_collider->add_endpoints();
                colliders.push_back(composite_collider);
                colliders.back()->id = colliders.size() - 1;
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = colliders.size() - 1;
                }
            }
            break;
        case MOBILE:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                if (dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]))
                {
                    console_log("SCRIPT ERROR: tried to make mobile collider from mobile collider");
                    return -1;
                }
                double translation_bounce = instruction->double_args[8];
                double rotation_bounce = instruction->double_args[9];
                if (translation_bounce < 1.0 || translation_bounce > 2.0 || rotation_bounce < 1.0 || rotation_bounce > 2.0)
                {
                    console_log("SCRIPT ERROR: invalid mobile collider motion bounce coefficient");
                    return -1;
                }
                MobileCollider* mobile_collider = new MobileCollider(colliders[instruction->int_args[0]], instruction->double_args[0], instruction->double_args[1], instruction->double_args[2], instruction->double_args[3], instruction->double_args[4], instruction->double_args[5], instruction->double_args[6], instruction->double_args[7]);
                mobile_collider->motion_bounce_coefficient = translation_bounce;
                mobile_collider->angular_motion_bounce_coefficient = rotation_bounce;
                colliders.push_back(mobile_collider);
                colliders.back()->id = colliders.size() - 1;
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    integer_variables[index] = colliders.size() - 1;
                }
            }
            break;
        case SOLID:
            if (table)
            {
                unsigned int layer = instruction->int_args[1];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                int collider_index = instruction->int_args[0];
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    if (instruction->int_args[0] == -1)
                    {
                        collider_index = colliders.size() - 1;
                    } else
                    {
                        console_log("SCRIPT ERROR: collider out of range");
                        return -1;
                    }
                }
                if (instruction->double_args[0] < 1.0 || instruction->double_args[0] > 2.0)
                {
                    console_log("SCRIPT ERROR: invalid bounce coefficient");
                    return -1;
                }
                if (instruction->double_args[1] < -1.0 || instruction->double_args[1] > 1.0)
                {
                    console_log("SCRIPT ERROR: invalid friction coefficient");
                    return -1;
                }
                Collider* collider = colliders[collider_index];
                collider->bounce_coefficient = instruction->double_args[0];
                collider->friction_coefficient = instruction->double_args[1];
                if (instruction->int_args[2] >= 0 && static_cast<unsigned int>(instruction->int_args[2]) < sounds.size())
                {
                    collider->impact_sound = sounds[instruction->int_args[2]];
                } else
                {
                    if (instruction->int_args[2] != -1)
                    {
                        console_log("SCRIPT ERROR: sound out of range");
                        return -1;
                    }
                }
                collider->solid = true;
                table->layers[layer]->add_collider(collider);
            } else
            {
                console_log("SCRIPT ERROR: no table declared");
                return -1;
            }
            break;
        case WOOD:
            if (table)
            {
                unsigned int layer = instruction->int_args[1];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                int collider_index = instruction->int_args[0];
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    if (instruction->int_args[0] == -1)
                    {
                        collider_index = colliders.size() - 1;
                    } else
                    {
                        console_log("SCRIPT ERROR: collider out of range");
                        return -1;
                    }
                }
                Collider* collider = colliders[collider_index];
                collider->bounce_coefficient = 1.25;
                collider->friction_coefficient = 0.44 * GLOBAL_FRICTION - 1.0;
                //TODO (Phase VII) default impact sounds
                collider->solid = true;
                table->layers[layer]->add_collider(collider);
            } else
            {
                console_log("SCRIPT ERROR: no table declared");
                return -1;
            }
            break;
        //TODO continue research on bounce coefficient for rubber posts and flippers; friction coefficients for all
        case METAL:
            if (table)
            {
                unsigned int layer = instruction->int_args[1];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                int collider_index = instruction->int_args[0];
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    if (instruction->int_args[0] == -1)
                    {
                        collider_index = colliders.size() - 1;
                    } else
                    {
                        console_log("SCRIPT ERROR: collider out of range");
                        return -1;
                    }
                }
                Collider* collider = colliders[collider_index];
                collider->bounce_coefficient = 1.15;
                collider->friction_coefficient = 0.22 * GLOBAL_FRICTION - 1.0;
                //TODO (Phase VII) default impact sounds
                collider->solid = true;
                table->layers[layer]->add_collider(collider);
            } else
            {
                console_log("SCRIPT ERROR: no table declared");
                return -1;
            }
            break;
        case PLASTIC:
            if (table)
            {
                unsigned int layer = instruction->int_args[1];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                int collider_index = instruction->int_args[0];
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    if (instruction->int_args[0] == -1)
                    {
                        collider_index = colliders.size() - 1;
                    } else
                    {
                        console_log("SCRIPT ERROR: collider out of range");
                        return -1;
                    }
                }
                Collider* collider = colliders[collider_index];
                collider->bounce_coefficient = 1.24;
                collider->friction_coefficient = 0.30 * GLOBAL_FRICTION - 1.0;
                //TODO (Phase VII) default impact sounds
                collider->solid = true;
                table->layers[layer]->add_collider(collider);
            } else
            {
                console_log("SCRIPT ERROR: no table declared");
                return -1;
            }
            break;
        case RUBBER:
            if (table)
            {
                unsigned int layer = instruction->int_args[1];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                int collider_index = instruction->int_args[0];
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    if (instruction->int_args[0] == -1)
                    {
                        collider_index = colliders.size() - 1;
                    } else
                    {
                        console_log("SCRIPT ERROR: collider out of range");
                        return -1;
                    }
                }
                Collider* collider = colliders[collider_index];
                //TODO more realistic test
                collider->bounce_coefficient = 1.47;
                collider->friction_coefficient = 1.92 * GLOBAL_FRICTION - 1.0;
                //TODO (Phase VII) default impact sounds
                collider->solid = true;
                table->layers[layer]->add_collider(collider);
            } else
            {
                console_log("SCRIPT ERROR: no table declared");
                return -1;
            }
            break;
        case RUBBER_BAND:
            if (table)
            {
                unsigned int layer = instruction->int_args[1];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                int collider_index = instruction->int_args[0];
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    if (instruction->int_args[0] == -1)
                    {
                        collider_index = colliders.size() - 1;
                    } else
                    {
                        console_log("SCRIPT ERROR: collider out of range");
                        return -1;
                    }
                }
                Collider* collider = colliders[collider_index];
                LineCollider* line_collider = dynamic_cast<LineCollider*>(collider);
                if (line_collider)
                {
                    line_collider->bounce_coefficient = 1.84;
                    //TODO (Phase VIII) allow band-type colliders to be declared generically
                    line_collider->band_bounce = true;
                    line_collider->outer_bounce_coefficient = 1.47; //TODO (Phase X) ensure that this value matches bounce coefficient of rubber posts
                    line_collider->friction_coefficient = 2.0 * GLOBAL_FRICTION - 1.0;
                    //TODO (Phase VII) default impact sounds
                    line_collider->solid = true;
                    table->layers[layer]->add_collider(line_collider);
                } else
                {
                    console_log("SCRIPT ERROR: non-line collider passed as line collider");
                    return -1;
                }
            } else
            {
                console_log("SCRIPT ERROR: no table declared");
                return -1;
            }
            break;
        case SWITCH:
            if (table)
            {
                unsigned int layer = instruction->int_args[1];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                int collider_index = instruction->int_args[0];
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    if (instruction->int_args[0] == -1)
                    {
                        collider_index = colliders.size() - 1;
                    } else
                    {
                        console_log("SCRIPT ERROR: collider out of range");
                        return -1;
                    }
                }
                Collider* collider = colliders[collider_index];
                collider->solid = false;
                collider->switch_id = next_switch_id;
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = next_switch_id;
                    ++next_switch_id;
                }
                table->layers[layer]->add_collider(collider);
            } else
            {
                console_log("SCRIPT ERROR: no table declared");
                return -1;
            }
            break;
        case GATE:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                Collider* collider = colliders[instruction->int_args[0]];
                collider->gate = instruction->int_args[1];
            }
            break;
        case MAGNET_POWER:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= magnets.size())
                {
                    console_log("SCRIPT ERROR: magnet out of range");
                    return -1;
                }
                Magnet* magnet = magnets[instruction->int_args[0]];
                magnet->magnetism = instruction->double_args[0];
            }
            break;
        case SYNC_MAGNET:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= magnets.size())
                {
                    console_log("SCRIPT ERROR: magnet out of range");
                    return -1;
                }
                Magnet* magnet = magnets[instruction->int_args[0]];
                if (instruction->int_args[1] < 0 || static_cast<unsigned int>(instruction->int_args[1]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[1]]);
                if (mobile_collider)
                {
                    mobile_collider->sync_magnet(magnet, instruction->double_args[0], instruction->double_args[1]);
                } else
                {
                    console_log("SCRIPT ERROR: attempted to sync magnet to non-mobile collider");
                    return -1;
                }
            }
            break;
        case LAYER_PORTAL:
            {
                //TODO (Phase IX) warn that layer portals must use collider ID, not switch ID
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                Collider* collider = colliders[instruction->int_args[0]];
                unsigned int layer = instruction->int_args[1];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                collider->portal = true;
                collider->destination_layer = layer;
            }
            break;
        case TANGIBLE:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                Collider* collider = colliders[instruction->int_args[0]];
                collider->solid = instruction->int_args[1];
            }
            break;
        case LIVE_CATCH:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                mobile_collider->live_catch_enabled = instruction->int_args[1];
            }
            break;
        case LIVE_CATCH_PROPERTIES:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                double elasticity = instruction->double_args[0];
                if (elasticity < 0.0)
                {
                    console_log("SCRIPT ERROR: invalid live catch elasticity");
                    return -1;
                }
                double damping = instruction->double_args[1];
                if (damping < 0.0 || damping > 1.0)
                {
                    console_log("SCRIPT ERROR: invalid live catch damping coefficient");
                    return -1;
                }
                mobile_collider->live_catch_elasticity = elasticity;
                mobile_collider->live_catch_damping_coefficient = damping;
            }
            break;
        case SET_SENSITIVITY:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                Collider* collider = colliders[instruction->int_args[0]];
                collider->sensitivity = instruction->double_args[0];
            }
            break;
        case SET_SPRITE_VISIBILITY:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= sprite_instances.size())
                {
                    console_log("SCRIPT ERROR: sprite out of range");
                    return -1;
                }
                SpriteInstance* sprite_instance = sprite_instances[instruction->int_args[0]];
                sprite_instance->visible = instruction->int_args[1];
            }
            break;
        case TRANSLATE_SPRITE:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= sprite_instances.size())
                {
                    console_log("SCRIPT ERROR: sprite out of range");
                    return -1;
                }
                SpriteInstance* sprite_instance = sprite_instances[instruction->int_args[0]];
                sprite_instance->position.x = instruction->double_args[0];
                sprite_instance->position.y = instruction->double_args[1];
            }
            break;
        case ROTATE_SPRITE:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= sprite_instances.size())
                {
                    console_log("SCRIPT ERROR: sprite out of range");
                    return -1;
                }
                SpriteInstance* sprite_instance = sprite_instances[instruction->int_args[0]];
                sprite_instance->angle = instruction->double_args[0];
            }
            break;
        case SET_BALL_SPRITE_ROTATION:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= balls.size())
                {
                    console_log("SCRIPT ERROR: ball out of range");
                    return -1;
                }
                Ball* ball = balls[instruction->int_args[0]];
                ball->sprite_rotation = instruction->int_args[1];
            }
            break;
        case SYNC_SPRITE:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= sprite_instances.size())
                {
                    console_log("SCRIPT ERROR: sprite out of range");
                    return -1;
                }
                SpriteInstance* sprite_instance = sprite_instances[instruction->int_args[0]];
                if (instruction->int_args[1] < 0 || static_cast<unsigned int>(instruction->int_args[1]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[1]]);
                if (mobile_collider)
                {
                    sprite_instance->sync(mobile_collider);
                    sprite_instance->synced_collider = mobile_collider;
                } else
                {
                    console_log("SCRIPT ERROR: attempted to sync sprite to non-mobile collider");
                    return -1;
                }
            }
            break;
        case SET_SPRITE_ANIMATION:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= sprite_instances.size())
                {
                    console_log("SCRIPT ERROR: sprite out of range");
                    return -1;
                }
                SpriteInstance* sprite_instance = sprite_instances[instruction->int_args[0]];
                sprite_instance->animation = instruction->int_args[1];
            }
            break;
        case SET_SPRITE_FRAME:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= sprite_instances.size())
                {
                    console_log("SCRIPT ERROR: sprite out of range");
                    return -1;
                }
                SpriteInstance* sprite_instance = sprite_instances[instruction->int_args[0]];
                sprite_instance->frame = instruction->int_args[1];
            }
            break;
        case LOOP_SPRITE_ANIMATION: //ID in, animation, start frame, frame count, frames per second, loop count
            {
                //TODO (Phase X) test (all arguments)
                //TODO (Phase VIII) validate arguments
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= sprite_instances.size())
                {
                    console_log("SCRIPT ERROR: sprite out of range");
                    return -1;
                }
                SpriteInstance* sprite_instance = sprite_instances[instruction->int_args[0]];
                sprite_instance->start_animation(instruction->int_args[1], instruction->int_args[3], instruction->double_args[0], instruction->int_args[4], 0, instruction->int_args[2], instruction->int_args[2]);
            }
            break;
        case LOOP_SYNCHRONOUS_SPRITE_ANIMATION: //ID in, animation, start frame, frame count, frames per second, sync offset
            {
                //TODO (Phase X) test (all arguments)
                //TODO (Phase VIII) validate arguments
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= sprite_instances.size())
                {
                    console_log("SCRIPT ERROR: sprite out of range");
                    return -1;
                }
                SpriteInstance* sprite_instance = sprite_instances[instruction->int_args[0]];
                sprite_instance->start_synchronous_animation(instruction->int_args[1], instruction->int_args[3], instruction->double_args[0], table->global_time + instruction->double_args[1], instruction->int_args[2]);
            }
            break;
        case NUDGE:
            {
                if (table)
                {
                    table->nudge(instruction->double_args[0], instruction->double_args[1], instruction->double_args[2]);
                } else
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
            }
            break;
        case HALT_NUDGE:
            {
                if (table)
                {
                    if (table->nudge_duration > 0.0)
                    {
                        table->nudge(-table->current_nudge.x, -table->current_nudge.y, 0.0);
                    }
                } else
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
            }
            break;
        case BUMPER_ACTION: //layer, collider ID in, added radius, speed
            {
                if (!table)
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
                unsigned int layer = instruction->int_args[0];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                int index = instruction->int_args[1];
                if (index < 0 || static_cast<unsigned int>(index) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                ArcCollider* collider = dynamic_cast<ArcCollider*>(colliders[index]);
                if (!collider)
                {
                    console_log("SCRIPT ERROR: non-arc collider passed as arc collider");
                }
                //TODO (Phase X) test multiple balls at once
                Point center = Point(collider->center.x, collider->center.y);
                double radius = collider->radius + instruction->double_args[0];
                double speed = instruction->double_args[1];
                for (unsigned int i = 0; i < table->layers[layer]->balls.size(); ++i)
                {
                    //TODO (Phase X) test to ensure that rapid contact between multiple bumpers will not break physics with excess speed
                    Ball* ball = table->layers[layer]->balls[i].first;
                    Point center_to_ball = Point(ball->position.x - center.x, ball->position.y - center.y);
                    double center_to_ball_square_magnitude = center_to_ball.square_magnitude();
                    if (center_to_ball_square_magnitude <= (radius + ball->radius) * (radius + ball->radius))
                    {
                        double center_to_ball_magnitude = center_to_ball.magnitude();
                        if (center_to_ball_magnitude)
                        {
                            Point center_to_ball_unit = Point(center_to_ball.x / center_to_ball_magnitude, center_to_ball.y / center_to_ball_magnitude);
                            ball->bounce(&center_to_ball_unit, collider->bounce_coefficient, collider->friction_coefficient);
                            ball->velocity.x += center_to_ball_unit.x * speed;
                            ball->velocity.y += center_to_ball_unit.y * speed;
                        }
                    }
                }
            }
            break;
        case SLINGSHOT_ACTION: //layer, collider ID in, kicker range, kicker speed, kicker bias
            {
                //TODO (Phase X) test multiple balls at once
                if (!table)
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
                unsigned int layer = instruction->int_args[0];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                int index = instruction->int_args[1];
                if (index < 0 || static_cast<unsigned int>(index) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                LineCollider* collider = dynamic_cast<LineCollider*>(colliders[index]);
                if (!collider)
                {
                    console_log("SCRIPT ERROR: non-line collider passed as line collider");
                }
                Point a = Point(collider->a.x, collider->a.y);
                Point b = Point(collider->b.x, collider->b.y);
                double range = instruction->double_args[0];
                if (range < 0.0)
                {
                    range *= -1.0;
                    std::swap(a.x, b.x);
                    std::swap(a.y, b.y);
                    //a_to_b.x = a.x - b.x;
                    //a_to_b.y = a.y - b.y;
                }
                Point a_to_b = Point(b.x - a.x, b.y - a.y);
                double a_to_b_magnitude = a_to_b.magnitude();
                Point a_to_b_unit = Point(a_to_b.x / a_to_b_magnitude, a_to_b.y / a_to_b_magnitude);
                double speed = instruction->double_args[1];
                double bias = instruction->double_args[2];
                if (bias < 0.0 || bias > 1.0)
                {
                    console_log("SCRIPT ERROR: invalid slingshot bias");
                    return -1;
                }
                for (unsigned int i = 0; i < table->layers[layer]->balls.size(); ++i)
                {
                    //TODO (Phase X) test negative range as opposite side slingshot
                    Ball* ball = table->layers[layer]->balls[i].first;
                    Point a_to_ball = Point(ball->position.x - a.x, ball->position.y - a.y);
                    double ball_bias = a_to_b_unit.dot_product(a_to_ball) / a_to_b_magnitude;
                    if (ball_bias > 0.0 && ball_bias < 1.0)
                    {
                        //TODO (Phase X) test bias system extensively
                        Point a_to_b_unit_normal = Point(-a_to_b_unit.y, a_to_b_unit.x);
                        double ball_distance = a_to_b_unit_normal.dot_product(a_to_ball) - ball->radius;
                        double bias_match = ball_bias > bias ? (1.0 - ball_bias) / (1.0 - bias) : ball_bias / bias;
                        double range_at_bias = range * bias_match;
                        if (ball_distance < range_at_bias)
                        {
                            double speed_at_bias = speed * std::sqrt(bias_match);
                            Point velocity_at_bias = Point(a_to_b_unit_normal.x * speed_at_bias, a_to_b_unit_normal.y * speed_at_bias);
                            //TODO (Phase X) ensure that the ball triggering the collider is not affected by this bounce
                            double bounce_coefficient = collider->bounce_coefficient;
                            if (collider->band_bounce)
                            {
                                double center_bias_match = ball_bias > 0.5 ? (1.0 - ball_bias) / 0.5 : ball_bias / 0.5;
                                double root_center_bias_match = std::sqrt(center_bias_match);
                                bounce_coefficient = root_center_bias_match * collider->bounce_coefficient + (1.0 - root_center_bias_match) * collider->outer_bounce_coefficient;
                            }
                            ball->bounce(&a_to_b_unit_normal, bounce_coefficient, collider->friction_coefficient, collider->impact_sound); //TODO apply bias here
                            ball->velocity.x += velocity_at_bias.x;
                            ball->velocity.y += velocity_at_bias.y;
                        }
                    }
                }
            }
            break;
        case KICKOUT_ACTION: //origin layer, x, y, destination layer, impulse x, impulse y
            {
                if (!table)
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
                unsigned int origin_layer = instruction->int_args[0];
                unsigned int destination_layer = instruction->int_args[1];
                if (origin_layer >= table->layers.size() || destination_layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                Point origin = Point(instruction->double_args[0], instruction->double_args[1]);
                Point velocity = Point(instruction->double_args[2], instruction->double_args[3]);
                for (unsigned int i = 0; i < table->layers[origin_layer]->balls.size(); ++i)
                {
                    Ball* ball = table->layers[origin_layer]->balls[i].first;
                    Point origin_to_ball = Point(ball->position.x - origin.x, ball->position.y - origin.y);
                    if (origin_to_ball.square_magnitude() < ball->radius * ball->radius)
                    {
                        //TODO (Phase X) test kicking onto other ball (including exactly on top, with and without speed)
                        table->layers[origin_layer]->pass_ball(ball, table->layers[destination_layer]);
                        ball->velocity.x = velocity.x;
                        ball->velocity.y = velocity.y;
                    }
                }
            }
            break;
        case SET_MOBILE_VELOCITY:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                mobile_collider->velocity.x = instruction->double_args[0];
                mobile_collider->velocity.y = instruction->double_args[1];
            }
            break;
        case SET_MOBILE_ACCELERATION:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                mobile_collider->acceleration.x = instruction->double_args[0];
                mobile_collider->acceleration.y = instruction->double_args[1];
            }
            break;
        case SET_MOBILE_ANGULAR_VELOCITY:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                mobile_collider->angular_velocity = instruction->double_args[0];
            }
            break;
        case SET_MOBILE_ANGULAR_ACCELERATION:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                mobile_collider->angular_acceleration = instruction->double_args[0];
            }
            break;
        case SET_MOBILE_SPRING_CONSTANT:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                mobile_collider->spring_constant = instruction->double_args[0];
            }
            break;
        case SET_MOBILE_ANGULAR_SPRING_CONSTANT:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                mobile_collider->angular_spring_constant = instruction->double_args[0];
            }
            break;
        case MOBILE_BOUNCE:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                mobile_collider->bounce_enabled = instruction->int_args[1];
            }
            break;
        case GET_MOBILE_POSITION:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                int x_index = instruction->int_args[1];
                if (x_index != -1)
                {
                    double_variables[x_index] = mobile_collider->origin.x + mobile_collider->offset.x;
                }
                int y_index = instruction->int_args[2];
                if (y_index != -1)
                {
                    double_variables[y_index] = mobile_collider->origin.y + mobile_collider->offset.y;
                }
            }
            break;
        case GET_MOBILE_VELOCITY:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                int x_index = instruction->int_args[1];
                if (x_index != -1)
                {
                    double_variables[x_index] = mobile_collider->velocity.x;
                }
                int y_index = instruction->int_args[2];
                if (y_index != -1)
                {
                    double_variables[y_index] = mobile_collider->velocity.y;
                }
            }
            break;
        case GET_MOBILE_ACCELERATION:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                int x_index = instruction->int_args[1];
                if (x_index != -1)
                {
                    double_variables[x_index] = mobile_collider->acceleration.x;
                }
                int y_index = instruction->int_args[2];
                if (y_index != -1)
                {
                    double_variables[y_index] = mobile_collider->acceleration.y;
                }
            }
            break;
        case GET_MOBILE_ANGLE:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    double_variables[index] = mobile_collider->angular_offset;
                }
            }
            break;
        case GET_MOBILE_ANGULAR_VELOCITY:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    double_variables[index] = mobile_collider->angular_velocity;
                }
            }
            break;
        case GET_MOBILE_ANGULAR_ACCELERATION:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= colliders.size())
                {
                    console_log("SCRIPT ERROR: collider out of range");
                    return -1;
                }
                MobileCollider* mobile_collider = dynamic_cast<MobileCollider*>(colliders[instruction->int_args[0]]);
                if (!mobile_collider)
                {
                    console_log("SCRIPT ERROR: non-mobile collider passed as mobile collider");
                    return -1;
                }
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    double_variables[index] = mobile_collider->angular_acceleration;
                }
            }
            break;
        case SET_BALL_POSITION:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= balls.size())
                {
                    console_log("SCRIPT ERROR: ball out of range");
                    return -1;
                }
                Ball* ball = balls[instruction->int_args[0]];
                ball->position.x = instruction->double_args[0];
                ball->position.y = instruction->double_args[1];
            }
            break;
        case SET_BALL_VELOCITY:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= balls.size())
                {
                    console_log("SCRIPT ERROR: ball out of range");
                    return -1;
                }
                Ball* ball = balls[instruction->int_args[0]];
                ball->velocity.x = instruction->double_args[0];
                ball->velocity.y = instruction->double_args[1];
            }
            break;
        case SET_BALL_LAYER:
            {
                if (table)
                {
                    if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= balls.size())
                    {
                        console_log("SCRIPT ERROR: ball out of range");
                        return -1;
                    }
                    unsigned int layer = instruction->int_args[1];
                    if (layer >= table->layers.size())
                    {
                        console_log("SCRIPT ERROR: layer out of range");
                        return -1;
                    }
                    Ball* ball = balls[instruction->int_args[0]];
                    Layer* origin = table->layers[ball->layer];
                    Layer* destination = table->layers[layer];
                    origin->pass_ball(ball, destination);
                } else
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
            }
            break;
        case SET_BALL_TANGIBILITY:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= balls.size())
                {
                    console_log("SCRIPT ERROR: ball out of range");
                    return -1;
                }
                Ball* ball = balls[instruction->int_args[0]];
                ball->tangible = instruction->int_args[1];
            }
            break;
        case GET_BALL_POSITION:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= balls.size())
                {
                    console_log("SCRIPT ERROR: ball out of range");
                    return -1;
                }
                Ball* ball = balls[instruction->int_args[0]];
                int x_index = instruction->int_args[1];
                int y_index = instruction->int_args[2];
                if (x_index != -1)
                {
                    double_variables[x_index] = ball->position.x;
                }
                if (y_index != -1)
                {
                    double_variables[y_index] = ball->position.y;
                }
            }
            break;
        case GET_BALL_VELOCITY:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= balls.size())
                {
                    console_log("SCRIPT ERROR: ball out of range");
                    return -1;
                }
                Ball* ball = balls[instruction->int_args[0]];
                int x_index = instruction->int_args[1];
                int y_index = instruction->int_args[2];
                if (x_index != -1)
                {
                    double_variables[x_index] = ball->velocity.x;
                }
                if (y_index != -1)
                {
                    double_variables[y_index] = ball->velocity.y;
                }
            }
            break;
        case GET_BALL_LAYER:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= balls.size())
                {
                    console_log("SCRIPT ERROR: ball out of range");
                    return -1;
                }
                Ball* ball = balls[instruction->int_args[0]];
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    integer_variables[index] = ball->layer;
                }
            }
            break;
        case GET_BALL_TANGIBILITY:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= balls.size())
                {
                    console_log("SCRIPT ERROR: ball out of range");
                    return -1;
                }
                Ball* ball = balls[instruction->int_args[0]];
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    integer_variables[index] = ball->tangible;
                }
            }
            break;
        case GET_DIGITAL_INPUT:
            {
                //TODO (Phase X) test
                Input input = static_cast<Input>(instruction->int_args[0]);
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    integer_variables[index] = input_handler->input_state[input].first;
                }
            }
            break;
        case GET_ANALOG_INPUT:
            {
                Input input = static_cast<Input>(instruction->int_args[0]);
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    double_variables[index] = input_handler->input_state[input].second;
                }
            }
            break;
        case DOT_MATRIX_DISPLAY:
            {
                if (instruction->int_args[5] < 0 || static_cast<unsigned int>(instruction->int_args[5]) >= sprites.size())
                {
                    console_log("SCRIPT ERROR: spritesheet out of range");
                    return -1;
                }
                Sprite* dot = sprites[instruction->int_args[5]];
                DotMatrixDisplay* dot_matrix_display = new DotMatrixDisplay(renderer, instruction->int_args[0], instruction->int_args[1], dot);
                dot_matrix_display->lit_r = instruction->int_args[2];
                dot_matrix_display->lit_g = instruction->int_args[3];
                dot_matrix_display->lit_b = instruction->int_args[4];
                dot_matrix_displays.push_back(dot_matrix_display);
                int index = instruction->int_args[6];
                if (index != -1)
                {
                    integer_variables[index] = dot_matrix_displays.size() - 1;
                }
            }
            break;
        case ALPHANUMERIC_DISPLAY:
            {
                char filename[1024];
                strcpy(filename, instruction->string_args[0].c_str());
                std::ifstream segment_file_test(filename);
                if (!segment_file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(filename);
                    console_log(error_message);
                    return -1;
                }
                char protocol_filename[1024];
                strcpy(protocol_filename, instruction->string_args[1].c_str());
                std::ifstream protocol_file_test(protocol_filename);
                if (!protocol_file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(protocol_filename);
                    console_log(error_message);
                    return -1;
                }
                AlphanumericDisplay* alphanumeric_display = new AlphanumericDisplay(renderer, filename, protocol_filename, instruction->int_args[5], instruction->int_args[6], instruction->int_args[1], instruction->int_args[0]);
                alphanumeric_display->lit_r = instruction->int_args[2];
                alphanumeric_display->lit_g = instruction->int_args[3];
                alphanumeric_display->lit_b = instruction->int_args[4];
                alphanumeric_displays.push_back(alphanumeric_display);
                int index = instruction->int_args[7];
                if (index != -1)
                {
                    integer_variables[index] = alphanumeric_displays.size() - 1;
                }
            }
            break;
        case SCORE_REEL:
            {
                char filename[1024];
                strcpy(filename, instruction->string_args[0].c_str());
                std::ifstream reel_file_test(filename);
                if (!reel_file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(filename);
                    console_log(error_message);
                }
                ScoreReel* score_reel = new ScoreReel(renderer, instruction->int_args[0], filename, instruction->int_args[1], instruction->int_args[2], instruction->int_args[3], instruction->double_args[0]);
                char roll_sound_filename[1024];
                if (instruction->string_args[1].size())
                {
                    strcpy(roll_sound_filename, instruction->string_args[1].c_str());
                    std::ifstream roll_sound_file_test(roll_sound_filename);
                    if (!roll_sound_file_test.good())
                    {
                        std::string error_message = "SCRIPT ERROR: can't read from ";
                        error_message.append(roll_sound_filename);
                        console_log(error_message);
                    }
                    score_reel->add_sound(roll_sound_filename);
                }
                score_reel->id = score_reels.size();
                score_reels.push_back(score_reel);
                bool invert_roll = instruction->int_args[4];
                score_reel->invert_roll = invert_roll;
                int index = instruction->int_args[5];
                if (index != -1)
                {
                    integer_variables[index] = score_reels.size() - 1;
                }
            }
            break;
        case DMD_SHOW_TEXT:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= dot_matrix_displays.size())
                {
                    console_log("SCRIPT ERROR: dot matrix display out of range");
                    return -1;
                }
                DotMatrixDisplay* dot_matrix_display = dot_matrix_displays[instruction->int_args[0]];
                dot_matrix_display->show_text(renderer, fonts[instruction->int_args[1]], instruction->string_args[0], instruction->int_args[2], instruction->int_args[3], static_cast<TextAlignment>(instruction->int_args[6]), instruction->int_args[4], instruction->int_args[5]);
            }
            break;
        case DMD_SHOW_NUMBER: //DMD ID in, number, font, x, y, max width, max height, alignment, decimal place, comma separation, trailing zero
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= dot_matrix_displays.size())
                {
                    console_log("SCRIPT ERROR: dot matrix display out of range");
                    return -1;
                }
                DotMatrixDisplay* dot_matrix_display = dot_matrix_displays[instruction->int_args[0]];
                if (instruction->int_args[2] < 0 || static_cast<unsigned int>(instruction->int_args[2]) >= fonts.size())
                {
                    console_log("SCRIPT ERROR: font out of range");
                    return -1;
                }
                dot_matrix_display->show_number(renderer, fonts[instruction->int_args[2]], instruction->int_args[1], instruction->int_args[3], instruction->int_args[4], static_cast<TextAlignment>(instruction->int_args[7]), instruction->int_args[5], instruction->int_args[6], instruction->int_args[9], instruction->int_args[10], instruction->int_args[8]);
            }
            break;
        case DMD_SHOW_SPRITE:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= dot_matrix_displays.size())
                {
                    console_log("SCRIPT ERROR: dot matrix display out of range");
                    return -1;
                }
                DotMatrixDisplay* dot_matrix_display = dot_matrix_displays[instruction->int_args[0]];
                if (instruction->int_args[1] < 0 || static_cast<unsigned int>(instruction->int_args[1]) >= sprites.size())
                {
                    console_log("SCRIPT ERROR: spritesheet out of range");
                    return -1;
                }
                Sprite* sprite = sprites[instruction->int_args[1]];
                int x = instruction->int_args[4];
                int y = instruction->int_args[5];
                TextAlignment alignment = static_cast<TextAlignment>(instruction->int_args[6]);
                switch (alignment)
                {
                case ALIGN_TOP_LEFT:
                    break;
                case ALIGN_TOP:
                    x -= sprite->width / 2;
                    break;
                case ALIGN_TOP_RIGHT:
                    x -= sprite->width;
                    break;
                case ALIGN_LEFT:
                    y -= sprite->height / 2;
                    break;
                case ALIGN_CENTER:
                    x -= sprite->width / 2;
                    y -= sprite->height / 2;
                    break;
                case ALIGN_RIGHT:
                    x -= sprite->width;
                    y -= sprite->height / 2;
                    break;
                case ALIGN_BOTTOM_LEFT:
                    y -= sprite->height;
                    break;
                case ALIGN_BOTTOM:
                    x -= sprite->width / 2;
                    y -= sprite->height;
                    break;
                case ALIGN_BOTTOM_RIGHT:
                    x -= sprite->width;
                    y -= sprite->height;
                    break;
                }
                dot_matrix_display->show_sprite(renderer, sprite, x, y, 0.0, instruction->int_args[2], instruction->int_args[3]);
            }
            break;
        case DMD_SHOW_LINE:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= dot_matrix_displays.size())
                {
                    console_log("SCRIPT ERROR: dot matrix display out of range");
                    return -1;
                }
                DotMatrixDisplay* dot_matrix_display = dot_matrix_displays[instruction->int_args[0]];
                dot_matrix_display->show_line(renderer, instruction->int_args[1], instruction->int_args[2], instruction->int_args[3], instruction->int_args[4], instruction->int_args[5], instruction->int_args[6]);
            }
            break;
        case DMD_SHOW_POINT:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= dot_matrix_displays.size())
                {
                    console_log("SCRIPT ERROR: dot matrix display out of range");
                    return -1;
                }
                DotMatrixDisplay* dot_matrix_display = dot_matrix_displays[instruction->int_args[0]];
                dot_matrix_display->show_point(renderer, instruction->int_args[1], instruction->int_args[2], instruction->int_args[3], instruction->int_args[4]);
            }
            break;
        case DMD_INVERT:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= dot_matrix_displays.size())
                {
                    console_log("SCRIPT ERROR: dot matrix display out of range");
                    return -1;
                }
                DotMatrixDisplay* dot_matrix_display = dot_matrix_displays[instruction->int_args[0]];
                dot_matrix_display->invert = instruction->int_args[1];
            }
            break;
        case DMD_CLEAR:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= dot_matrix_displays.size())
                {
                    console_log("SCRIPT ERROR: dot matrix display out of range");
                    return -1;
                }
                DotMatrixDisplay* dot_matrix_display = dot_matrix_displays[instruction->int_args[0]];
                dot_matrix_display->clear_display(renderer, instruction->int_args[1]);
            }
            break;
        case AND_SHOW_TEXT:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= alphanumeric_displays.size())
                {
                    console_log("SCRIPT ERROR: alphanumeric display out of range");
                    return -1;
                }
                AlphanumericDisplay* alphanumeric_display = alphanumeric_displays[instruction->int_args[0]];
                alphanumeric_display->show_text(instruction->string_args[0], instruction->int_args[1], static_cast<TextAlignment>(instruction->int_args[2]), 0);
            }
            break;
        case AND_SHOW_NUMBER: //AND ID in, number, anchor character, alignment, decimal place, comma separation, trailing zero
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= alphanumeric_displays.size())
                {
                    console_log("SCRIPT ERROR: alphanumeric display out of range");
                    return -1;
                }
                AlphanumericDisplay* alphanumeric_display = alphanumeric_displays[instruction->int_args[0]];
                alphanumeric_display->show_number(instruction->int_args[1], instruction->int_args[2], static_cast<TextAlignment>(instruction->int_args[3]), 0, instruction->int_args[5], instruction->int_args[6], instruction->int_args[4]);
            }
            break;
        case AND_SET_SEGMENT:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= alphanumeric_displays.size())
                {
                    console_log("SCRIPT ERROR: alphanumeric display out of range");
                    return -1;
                }
                AlphanumericDisplay* alphanumeric_display = alphanumeric_displays[instruction->int_args[0]];
                alphanumeric_display->set_segment(instruction->int_args[1], instruction->int_args[2], instruction->int_args[3]);
            }
            break;
        case AND_CLEAR:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= alphanumeric_displays.size())
                {
                    console_log("SCRIPT ERROR: alphanumeric display out of range");
                    return -1;
                }
                AlphanumericDisplay* alphanumeric_display = alphanumeric_displays[instruction->int_args[0]];
                alphanumeric_display->clear_display(instruction->int_args[1]);
            }
            break;
        case REEL_SHOW_NUMBER:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= score_reels.size())
                {
                    console_log("SCRIPT ERROR: score reel out of range");
                    return -1;
                }
                ScoreReel* score_reel = score_reels[instruction->int_args[0]];
                score_reel->show_number(instruction->int_args[1]);
            }
            break;
        case REEL_RESET:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= score_reels.size())
                {
                    console_log("SCRIPT ERROR: score reel out of range");
                    return -1;
                }
                ScoreReel* score_reel = score_reels[instruction->int_args[0]];
                score_reel->show_number(0);
            }
            break;
        case REEL_ADD_SHADING:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= score_reels.size())
                {
                    console_log("SCRIPT ERROR: score reel out of range");
                    return -1;
                }
                ScoreReel* score_reel = score_reels[instruction->int_args[0]];
                std::string filename = instruction->string_args[0];
                std::ifstream file_test(filename);
                if (!file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(filename);
                    console_log(error_message);
                    return -1;
                    //TODO check all file test errors for missing return statement
                }
                score_reel->add_shading_texture(renderer, filename.c_str());
            }
            break;
        case REEL_LIGHT:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= score_reels.size())
                {
                    console_log("SCRIPT ERROR: score reel out of range");
                    return -1;
                }
                ScoreReel* score_reel = score_reels[instruction->int_args[0]];
                score_reel->lit = instruction->int_args[1];
            }
            break;
        case LOAD_FONT_SIZE:
            {
                char filename[1024];
                strcpy(filename, instruction->string_args[0].c_str());
                std::ifstream file_test(filename);
                if (!file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(filename);
                    console_log(error_message);
                    return -1;
                }
                char kerning_data_filename[1024];
                strcpy(kerning_data_filename, instruction->string_args[1].c_str());
                std::ifstream kerning_data_file_test(kerning_data_filename);
                if (!kerning_data_file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(kerning_data_filename);
                    console_log(error_message);
                    return -1;
                }
                FontSize* font_size = new FontSize(renderer, filename, kerning_data_filename);
                font_sizes.push_back(font_size);
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = font_sizes.size() - 1;
                }
            }
            break;
        case LOAD_MONOSPACED_FONT_SIZE:
            {
                char filename[1024];
                strcpy(filename, instruction->string_args[0].c_str());
                std::ifstream file_test(filename);
                if (!file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(filename);
                    console_log(error_message);
                    return -1;
                }
                FontSize* font_size = new FontSize(renderer, filename, instruction->int_args[0], instruction->int_args[1]);
                font_sizes.push_back(font_size);
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = font_sizes.size() - 1;
                }
            }
            break;
        case CREATE_FONT:
            {
                Font* font = new Font();
                fonts.push_back(font);
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = fonts.size() - 1;
                }
            }
            break;
        case ADD_FONT_SIZE:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= fonts.size())
                {
                    console_log("SCRIPT ERROR: font out of range");
                    return -1;
                }
                Font* font = fonts[instruction->int_args[0]];
                if (instruction->int_args[1] < 0 || static_cast<unsigned int>(instruction->int_args[1]) >= font_sizes.size())
                {
                    console_log("SCRIPT ERROR: font size out of range");
                    return -1;
                }
                FontSize* font_size = font_sizes[instruction->int_args[1]];
                font->add_size(font_size);
            }
            break;
        case BACKGLASS:
            {
                Backglass* backglass = new Backglass();
                backglass->top_left.x = instruction->double_args[0];
                backglass->top_left.y = instruction->double_args[1];
                backglass->dimensions.x = instruction->double_args[2];
                backglass->dimensions.y = instruction->double_args[3];
                backglasses.push_back(backglass);
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = backglasses.size() - 1;
                }
            }
            break;
        case ADD_DOT_MATRIX_DISPLAY:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= backglasses.size())
                {
                    console_log("SCRIPT ERROR: backglass out of range");
                    return -1;
                }
                Backglass* backglass = backglasses[instruction->int_args[0]];
                if (instruction->int_args[1] < 0 || static_cast<unsigned int>(instruction->int_args[1]) >= dot_matrix_displays.size())
                {
                    console_log("SCRIPT ERROR: dot matrix display out of range");
                    return -1;
                }
                DotMatrixDisplay* dot_matrix_display = dot_matrix_displays[instruction->int_args[1]];
                backglass->add_element(dot_matrix_display, instruction->double_args[0], instruction->double_args[1]);
            }
            break;
        case ADD_ALPHANUMERIC_DISPLAY:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= backglasses.size())
                {
                    console_log("SCRIPT ERROR: backglass out of range");
                    return -1;
                }
                Backglass* backglass = backglasses[instruction->int_args[0]];
                if (instruction->int_args[1] < 0 || static_cast<unsigned int>(instruction->int_args[1]) >= alphanumeric_displays.size())
                {
                    console_log("SCRIPT ERROR: alphanumeric display out of range");
                    return -1;
                }
                AlphanumericDisplay* alphanumeric_display = alphanumeric_displays[instruction->int_args[1]];
                backglass->add_element(alphanumeric_display, instruction->double_args[0], instruction->double_args[1]);
            }
            break;
        case ADD_SCORE_REEL:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= backglasses.size())
                {
                    console_log("SCRIPT ERROR: backglass out of range");
                    return -1;
                }
                Backglass* backglass = backglasses[instruction->int_args[0]];
                if (instruction->int_args[1] < 0 || static_cast<unsigned int>(instruction->int_args[1]) >= score_reels.size())
                {
                    console_log("SCRIPT ERROR: score reel out of range");
                    return -1;
                }
                ScoreReel* score_reel = score_reels[instruction->int_args[1]];
                backglass->add_element(score_reel, instruction->double_args[0], instruction->double_args[1]);
            }
            break;
        case DESIGNATE_BACKGLASS:
            {
                if (table)
                {
                    if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= backglasses.size())
                    {
                        console_log("SCRIPT ERROR: backglass out of range");
                        return -1;
                    }
                    Backglass* backglass = backglasses[instruction->int_args[0]];
                    table->backglass = backglass;
                } else
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
            }
            break;
        case EMBED_BACKGLASS:
            {
                if (table)
                {
                    unsigned int layer = instruction->int_args[1];
                    if (layer >= table->layers.size())
                    {
                        console_log("SCRIPT ERROR: layer out of range");
                        return -1;
                    }
                    if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= backglasses.size())
                    {
                        console_log("SCRIPT ERROR: backglass out of range");
                        return -1;
                    }
                    Backglass* backglass = backglasses[instruction->int_args[0]];
                    table->layers[layer]->add_backglass(backglass);
                } else
                {
                    console_log("SCRIPT ERROR: no table declared");
                    return -1;
                }
            }
            break;
        case PLAY_SOUND:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= sounds.size())
                {
                    console_log("SCRIPT ERROR: sound out of range");
                    return -1;
                }
                Sound* sound = sounds[instruction->int_args[0]];
                sound->play(instruction->int_args[1]);
            }
            break;
        case STOP_SOUND:
            {
                Mix_HaltGroup(TRACK_TABLE);
            }
            break;
        case PLAY_MUSIC:
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= songs.size())
                {
                    console_log("SCRIPT ERROR: song out of range");
                    return -1;
                }
                Mix_Music* song = songs[instruction->int_args[0]];
                Mix_VolumeMusic(instruction->int_args[1]);
                Mix_PlayMusic(song, -1);
                song_playing = song;
            }
            break;
        case TRANSITION_MUSIC:
            {
                //TODO (Phase X) test
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= songs.size())
                {
                    console_log("SCRIPT ERROR: song out of range");
                    return -1;
                }
                Mix_Music* song = songs[instruction->int_args[0]];
                Mix_VolumeMusic(instruction->int_args[1]);
                double position = 0.0;
                if (song_playing)
                {
                    Mix_GetMusicPosition(song_playing);
                }
                Mix_PlayMusic(song, -1);
                Mix_SetMusicPosition(position);
                song_playing = song;
            }
            break;
        case STOP_MUSIC:
            {
                Mix_HaltMusic();
            }
            break;
        case INTEGER:
            {
                std::string name = instruction->string_args[0];
                if (syntax->enum_table.find(name) != syntax->enum_table.end() || name == "NULL")
                {
                    console_log("SCRIPT ERROR: variable name reserved");
                    return -1;
                }
                if (variable_index_table.find(name) == variable_index_table.end())
                {
                    integer_variables.push_back(instruction->int_args[0]);
                    variable_index_table[name] = std::pair<DataType, int>(TYPE_INTEGER_VARIABLE, integer_variables.size() - 1);
                } else
                {
                    integer_variables[variable_index_table[name].second] = instruction->int_args[0];
                }
            }
            break;
        case DOUBLE:
            {
                std::string name = instruction->string_args[0];
                if (syntax->enum_table.find(name) != syntax->enum_table.end() || name == "NULL")
                {
                    console_log("SCRIPT ERROR: variable name reserved");
                    return -1;
                }
                if (variable_index_table.find(name) == variable_index_table.end())
                {
                    double_variables.push_back(instruction->double_args[0]);
                    variable_index_table[name] = std::pair<DataType, int>(TYPE_DOUBLE_VARIABLE, double_variables.size() - 1);
                } else
                {
                    double_variables[variable_index_table[name].second] = instruction->double_args[0];
                }
            }
            break;
        case STRING:
            {
                std::string name = instruction->string_args[0];
                if (syntax->enum_table.find(name) != syntax->enum_table.end() || name == "NULL")
                {
                    console_log("SCRIPT ERROR: variable name reserved");
                    return -1;
                }
                if (variable_index_table.find(name) == variable_index_table.end())
                {
                    string_variables.push_back(instruction->string_args[1]);
                    variable_index_table[name] = std::pair<DataType, int>(TYPE_STRING_VARIABLE, string_variables.size() - 1);
                } else
                {
                    string_variables[variable_index_table[name].second] = instruction->string_args[1];
                }
            }
            break;
        case ARRAY:
            {
                std::string name = instruction->string_args[0];
                if (syntax->enum_table.find(name) != syntax->enum_table.end() || name == "NULL")
                {
                    console_log("SCRIPT ERROR: array name reserved");
                    return -1;
                }
                int n = instruction->int_args[0];
                int value = instruction->int_args[1];
                if (variable_index_table.find(name) == variable_index_table.end())
                {
                    int* start = new int[n];
                    for (int i = 0; i < n; ++i)
                    {
                        start[i] = value;
                    }
                    arrays.push_back(std::pair<int, int*>(n, start));
                    variable_index_table[name] = std::pair<DataType, int>(TYPE_ARRAY, arrays.size() - 1);
                } else
                {
                    delete arrays[variable_index_table[name].second].second;
                    int* start = new int[n];
                    for (int i = 0; i < n; ++i)
                    {
                        start[i] = value;
                    }
                    arrays[variable_index_table[name].second] = std::pair<int, int*>(n, start);
                }
            }
            break;
        case SET_INTEGER:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[1];
                }
            }
            break;
        case SET_DOUBLE:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    double_variables[index] = instruction->double_args[0];
                }
            }
            break;
        case SET_STRING:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    string_variables[index] = instruction->string_args[0];
                }
            }
            break;
            //TODO (Phase VIII) abort when bounds check fails
        case SET_ARRAY_ELEMENT:
            {
                int index = instruction->int_args[0];
                int inner_index = instruction->int_args[1];
                int value = instruction->int_args[2];
                if (index != -1)
                {
                    int n = arrays[index].first;
                    if (inner_index >= 0 && inner_index < n)
                    {
                        int* start = arrays[index].second;
                        start[inner_index] = value;
                    } else
                    {
                        console_log("SCRIPT ERROR: array element out of range");
                        return -1;
                    }
                }
            }
            break;
        case GET_ARRAY_ELEMENT:
            {
                int index = instruction->int_args[0];
                int inner_index = instruction->int_args[1];
                int index_out = instruction->int_args[2];
                if (index != -1)
                {
                    int n = arrays[index].first;
                    if (inner_index >= 0 && inner_index < n)
                    {
                        if (index_out != -1)
                        {
                            int* start = arrays[index].second;
                            int value = start[inner_index];
                            integer_variables[index_out] = value;
                        }
                    } else
                    {
                        console_log("SCRIPT ERROR: array element out of range");
                        return -1;
                    }
                }
            }
            break;
        case STRING_LENGTH:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = instruction->string_args[0].size();
                }
            }
            break;
        case STRING_CHARACTER:
            {
                int char_index = instruction->int_args[0];
                if (char_index < 0 || static_cast<unsigned int>(char_index) >= instruction->string_args[0].size())
                {
                    console_log("SCRIPT ERROR: character out of range");
                    return -1;
                }
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    integer_variables[index] = instruction->string_args[0][char_index];
                }
            }
            break;
        case STRING_APPEND:
            {
                //TODO (Phase X) test to confirm that STRING_VARIABLE arguments are processed as integer indices
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    string_variables[index].push_back(instruction->int_args[0]);
                }
            }
            break;
        case STRING_CONCATENATE:
            {
                //TODO (Phase X) test
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    string_variables[index].append(instruction->string_args[0]);
                }
            }
            break;
        case RANDOM:
            {
                int range = instruction->int_args[0];
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    integer_variables[index] = std::rand() % range;
                }
            }
            break;
        case ADD:
            {
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] + instruction->int_args[1];
                }
            }
            break;
        case SUBTRACT:
            {
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] - instruction->int_args[1];
                }
            }
            break;
        case MULTIPLY:
            {
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] * instruction->int_args[1];
                }
            }
            break;
        case DIVIDE:
            {
                if (instruction->int_args[1] == 0)
                {
                    console_log("SCRIPT ERROR: tried to divide by zero");
                    return -1;
                }
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] / instruction->int_args[1];
                }
            }
            break;
        case MODULO:
            {
                if (instruction->int_args[1] == 0)
                {
                    console_log("SCRIPT ERROR: tried to divide by zero");
                    return -1;
                }
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] % instruction->int_args[1];
                }
            }
            break;
        case D_ADD:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    double_variables[index] = instruction->double_args[0] + instruction->double_args[1];
                }
            }
            break;
        case D_SUBTRACT:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    double_variables[index] = instruction->double_args[0] - instruction->double_args[1];
                }
            }
            break;
        case D_MULTIPLY:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    double_variables[index] = instruction->double_args[0] * instruction->double_args[1];
                }
            }
            break;
        case D_DIVIDE:
            {
                if (instruction->double_args[1] == 0.0)
                {
                    console_log("SCRIPT ERROR: tried to divide by zero");
                    return -1;
                }
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    double_variables[index] = instruction->double_args[0] / instruction->double_args[1];
                }
            }
            break;
        case AND:
            {
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] && instruction->int_args[1];
                }
            }
            break;
        case OR:
            {
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] || instruction->int_args[1];
                }
            }
            break;
        case XOR:
            {
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = (instruction->int_args[0] != 0) ^ (instruction->int_args[1] != 0);
                }
            }
            break;
        case NOT:
            {
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    integer_variables[index] = !instruction->int_args[0];
                }
            }
            break;
        case BITWISE_AND:
            {
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] & instruction->int_args[1];
                }
            }
            break;
        case BITWISE_OR:
            {
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] | instruction->int_args[1];
                }
            }
            break;
        case BITWISE_XOR:
            {
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] ^ instruction->int_args[1];
                }
            }
            break;
        case BITWISE_NOT:
            {
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    integer_variables[index] = ~instruction->int_args[0];
                }
            }
            break;
        case GREATER:
            {
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] > instruction->int_args[1];
                }
            }
            break;
        case LESS:
            {
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] < instruction->int_args[1];
                }
            }
            break;
        case EQUAL:
            {
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = instruction->int_args[0] == instruction->int_args[1];
                }
            }
            break;
        case D_GREATER:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = instruction->double_args[0] > instruction->double_args[1];
                }
            }
            break;
        case D_LESS:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = instruction->double_args[0] < instruction->double_args[1];
                }
            }
            break;
        case D_EQUAL:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = instruction->double_args[0] == instruction->double_args[1];
                }
            }
            break;
        case SHIFT_ARRAY:
            {
                //TODO (Phase X) test thoroughly (positive and negative shifts, small and large shifts)
                int index = instruction->int_args[0];
                int shift_length = instruction->int_args[1];
                int array_length = arrays[index].first;
                int* to_shift = arrays[index].second;
                int* temp_array = new int[array_length];
                for (int i = 0; i < array_length; ++i)
                {
                    temp_array[(i + shift_length + array_length) % array_length] = to_shift[i];
                }
                for (int i = 0; i < array_length; ++i)
                {
                    to_shift[i] = temp_array[i];
                }
                delete[] temp_array;
            }
            break;
        case CLEAR_ARRAY:
            {
                //TODO bounds checking
                int index = instruction->int_args[0];
                for (int i = 0; i < arrays[index].first; ++i)
                {
                    arrays[index].second[i] = instruction->int_args[1];
                }
            }
            break;
        case INTEGER_TO_DOUBLE:
            {
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    double_variables[index] = 1.0 * instruction->int_args[0];
                }
            }
            break;
        case DOUBLE_TO_INTEGER:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = std::floor(instruction->double_args[0]);
                }
            }
            break;
        case INTEGER_TO_STRING:
            {
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    string_variables[index] = std::to_string(instruction->int_args[0]);
                }
            }
            break;
        case DOUBLE_TO_STRING:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    string_variables[index] = std::to_string(instruction->double_args[0]);
                }
            }
            break;
            //TODO (Phase X) test timer commands
            //TODO bounds checking on timer commands
        case TIMER:
            {
                std::string name = instruction->string_args[0];
                if (syntax->enum_table.find(name) != syntax->enum_table.end() || name == "NULL")
                {
                    console_log("SCRIPT ERROR: variable name reserved");
                    return -1;
                }
                if (variable_index_table.find(name) == variable_index_table.end())
                {
                    integer_variables.push_back(0);
                    variable_index_table[name] = std::pair<DataType, int>(TYPE_INTEGER_VARIABLE, integer_variables.size() - 1);
                } else
                {
                    integer_variables[variable_index_table[name].second] = 0;
                }
                timers[integer_variables.size() - 1] = false;
            }
            break;
        case START_TIMER:
            {
                //TODO consider making START_TIMER 0.0 instantly queue TIMER_DONE
                int index = instruction->int_args[0];
                int time = instruction->double_args[0] * PHYSICS_TICKS_PER_FRAME * FRAMES_PER_SECOND;
                if (index != -1)
                {
                    integer_variables[index] = time;
                    timers[index] = true;
                }
            }
            break;
        case PAUSE_TIMER:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    timers[index] = false;
                }
            }
            break;
        case RESUME_TIMER:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    timers[index] = true;
                }
            }
            break;
        case TIMER_TO_DOUBLE:
            {
                int timer_index = instruction->int_args[0];
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    double_variables[index] = integer_variables[timer_index] * 1.0 / (PHYSICS_TICKS_PER_FRAME * FRAMES_PER_SECOND);
                }
            }
            break;
        case IF:
            //Only for flow control
            break;
        case ELSE:
            //Only for flow control
            break;
        case ENDIF:
            //Only for flow control
            break;
        case WHILE:
            //Only for flow control
            break;
        case WEND:
            //Only for flow control
            break;
        case BREAK:
            //Only for flow control
            break;
        case GOSUB:
            //Only for flow control
            break;
        case RETURN:
            //Only for flow control
            break;
        case LABEL:
            //Handled at parse time
            break;
        case LOAD_SPRITESHEET:
            {
                char filename[1024];
                strcpy(filename, instruction->string_args[0].c_str());
                std::ifstream file_test(filename);
                if (!file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(filename);
                    console_log(error_message);
                    return -1;
                }
                Sprite* sprite = new Sprite(renderer, filename, instruction->double_args[0], instruction->double_args[1], instruction->int_args[0], instruction->int_args[1]);
                sprites.push_back(sprite);
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = sprites.size() - 1;
                }
            }
            break;
        case LOAD_SOUND:
            {
                char filename[1024];
                strcpy(filename, instruction->string_args[0].c_str());
                std::ifstream file_test(filename);
                if (!file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(filename);
                    console_log(error_message);
                    return -1;
                }
                Sound* sound = new Sound(filename, static_cast<AudioTrack>(instruction->int_args[0]));
                sounds.push_back(sound);
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    integer_variables[index] = sounds.size() - 1;
                }
            }
            break;
        case LOAD_MUSIC:
            {
                //TODO (Phase X) test
                char filename[1024];
                strcpy(filename, instruction->string_args[0].c_str());
                std::ifstream file_test(filename);
                if (!file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(filename);
                    console_log(error_message);
                    return -1;
                }
                Mix_Music* song = Mix_LoadMUS(filename);
                songs.push_back(song);
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = songs.size() - 1;
                }
            }
            break;
        case LOAD_RECORD: //filename, entry count, default name, default score, ID out
            {
                char filename[1024];
                strcpy(filename, instruction->string_args[0].c_str());
                /*
                std::ifstream file_test(filename);
                if (!file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(filename);
                    console_log(error_message);
                    return -1;
                }
                */
                records.push_back(new Record(filename, instruction->int_args[0], instruction->string_args[1], instruction->int_args[1]));
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = records.size() - 1;
                }
            }
            break;
        case LOAD_CONFIG: //filename
            {
                load_config(instruction->string_args[0].c_str(), this);
            }
            break;
        case LOAD_CONTROLS: //filename
            {
                char filename[1024];
                strcpy(filename, instruction->string_args[0].c_str());
                std::ifstream file_test(filename);
                if (!file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(filename);
                    console_log(error_message);
                    return -1;
                }
                input_handler->load_controls(filename);
                CONTROLS_PATH = filename;
            }
            break;
        case SAVE_RECORD: //ID in, filename, name, score
            {
                //TODO debug (currently appears to overwrite higher scores with lower scores
                //TODO guard against overwriting non-record files
                char filename[1024];
                strcpy(filename, instruction->string_args[0].c_str());
                /*
                std::ifstream file_test(filename);
                if (!file_test.good())
                {
                    std::string error_message = "SCRIPT ERROR: can't read from ";
                    error_message.append(filename);
                    console_log(error_message);
                    return -1;
                }
                */
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= records.size())
                {
                    console_log("SCRIPT ERROR: record out of range");
                    return -1;
                }
                Record* record = records[instruction->int_args[0]];
                std::string name = instruction->string_args[1];
                if (name == "")
                {
                    name = "_";
                }
                int score = instruction->int_args[1];
                if (record->qualifies(score))
                {
                    record->add_record(name, score);
                    record->save(filename);
                }
            }
            break;
        case SAVE_CONFIG: //filename
            {
                save_config(instruction->string_args[0].c_str(), this);
            }
            break;
        case SAVE_CONTROLS: //filename
            {
                char filename[1024];
                strcpy(filename, instruction->string_args[0].c_str());
                std::ifstream file_test(filename);
                input_handler->save_controls(filename);
            }
            break;
        case GET_RECORD: //ID in, rank, name out, score out
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= records.size())
                {
                    console_log("SCRIPT ERROR: record out of range");
                    return -1;
                }
                Record* record = records[instruction->int_args[0]];
                if (instruction->int_args[1] - 1 || static_cast<unsigned int>(instruction->int_args[1] - 1) >= record->scores.size())
                {
                    console_log("SCRIPT ERROR: record rank out of range");
                    return -1;
                }
                int name_index = instruction->int_args[2];
                int score_index = instruction->int_args[3];
                if (name_index != -1)
                {
                    string_variables[name_index] = record->scores[instruction->int_args[1] - 1].first;
                }
                if (score_index != -1)
                {
                    integer_variables[score_index] = record->scores[instruction->int_args[1] - 1].second;
                }
            }
            break;
        case QUALIFIES: //ID in, score, rank out (0 if not qualified)
            {
                if (instruction->int_args[0] < 0 || static_cast<unsigned int>(instruction->int_args[0]) >= records.size())
                {
                    console_log("SCRIPT ERROR: record out of range");
                    return -1;
                }
                Record* record = records[instruction->int_args[0]];
                int score = instruction->int_args[1];
                int index = instruction->int_args[2];
                if (index != -1)
                {
                    integer_variables[index] = record->qualifies(score);
                }
            }
            break;
        case DEFAULT_CONFIG:
            {
                default_config(this);
            }
            break;
        case DEFAULT_CONTROLS:
            {
                input_handler->default_controls();
            }
            break;
        case BIND_KEY: //key, input, negative edge
            {
                input_handler->bind_key(instruction->int_args[0], static_cast<Input>(instruction->int_args[1]), instruction->int_args[2], 1.0);
            }
            break;
        case BIND_BUTTON: //button, input, negative edge
            {
                input_handler->bind_button(static_cast<SDL_GameControllerButton>(instruction->int_args[0]), static_cast<Input>(instruction->int_args[1]), instruction->int_args[2], 1.0);
            }
            break;
        case BIND_AXIS: //axis, input, invert, dead zone, sensitivity
            {
                input_handler->bind_axis(static_cast<SDL_GameControllerAxis>(instruction->int_args[0]), static_cast<Input>(instruction->int_args[1]), instruction->int_args[2], instruction->double_args[0], instruction->double_args[1]);
            }
            break;
        case UNBIND_KEY: //key
            {
                auto bind_list = input_handler->keyboard_binds[instruction->int_args[0]];
                while (bind_list.size())
                {
                    bind_list.pop_back();
                }
            }
            break;
        case UNBIND_BUTTON: //button
            {
                auto bind_list = input_handler->gamepad_button_binds[static_cast<SDL_GameControllerButton>(instruction->int_args[0])];
                while (bind_list.size())
                {
                    bind_list.pop_back();
                }
            }
            break;
        case UNBIND_AXIS: //axis
            {
                auto bind_list = input_handler->gamepad_axis_binds[static_cast<SDL_GameControllerAxis>(instruction->int_args[0])];
                while (bind_list.size())
                {
                    bind_list.pop_back();
                }
            }
            break;
        case UNBIND_ALL:
            {
                input_handler->keyboard_binds.clear();
                input_handler->gamepad_button_binds.clear();
                input_handler->gamepad_axis_binds.clear();
            }
            break;
        case DISPLAY_MODE:
            {
                int mode = instruction->int_args[0];
                FULLSCREEN = mode;
                SDL_SetWindowFullscreen(window, mode ? SDL_WINDOW_FULLSCREEN_DESKTOP : 0);
                if (!FULLSCREEN)
                {
                    mode = instruction->int_args[1];
                    BORDERLESS = mode;
                    SDL_SetWindowBordered(window, !mode ? SDL_TRUE : SDL_FALSE);
                }
            }
            break;
        case RESOLUTION: //width, height
            {
                if (instruction->int_args[0] < 0 || instruction->int_args[1] < 0)
                {
                    console_log("SCRIPT ERROR: invalid resolution");
                    return -1;
                }
                SDL_SetWindowSize(window, instruction->int_args[0], instruction->int_args[1]);
                SDL_SetWindowPosition(window, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED);
            }
            break;
        case FRAMERATE:
            {
                int framerate = instruction->int_args[0];
                if (!framerate || 1440 % framerate)
                {
                    console_log("SCRIPT ERROR: invalid framerate (must be a factor of 1440)");
                } else
                {
                    FRAMES_PER_SECOND = framerate;
                    MS_PER_FRAME = 1000.0 / framerate;
                    PHYSICS_TICKS_PER_FRAME = 24 * 60 / framerate;
                }
            }
            break;
        case FLASHING_LIGHTS:
            {
                DISABLE_FLASHING = !instruction->int_args[0];
            }
            break;
        case FLASHING_LIGHTS_ENABLED:
            {
                int index = instruction->int_args[0];
                if (index != -1)
                {
                    integer_variables[index] = !DISABLE_FLASHING;
                }
            }
            break;
        case MUTE:
            {
                mute = instruction->int_args[0];
                if (mute)
                {
                    Mix_MasterVolume(0);
                    Mix_VolumeMusic(0);
                    Mix_Volume(-1, 0);
                } else
                {
                    Mix_MasterVolume(volume);
                    Mix_VolumeMusic(music_volume);
                    Mix_Volume(-1, sound_volume);
                }
            }
            break;
        case MASTER_VOLUME:
            {
                volume = instruction->int_args[0];
                if (!mute)
                {
                    Mix_MasterVolume(volume);
                }
            }
            break;
        case MUSIC_VOLUME:
            {
                music_volume = instruction->int_args[0];
                if (!mute)
                {
                    Mix_VolumeMusic(music_volume);
                }
            }
            break;
        case SOUND_VOLUME:
            {
                sound_volume = instruction->int_args[0];
                if (!mute)
                {
                    Mix_Volume(-1, sound_volume);
                }
            }
            break;
        case SCALE_UI:
            {
                if (instruction->double_args[0] < 0.0)
                {
                    console_log("SCRIPT ERROR: tried to set scale to a negative number");
                }
                UI_SCALE = instruction->double_args[0];
            }
            break;
        case SCALE_HUD:
            {
                if (instruction->double_args[0] < 0.0)
                {
                    console_log("SCRIPT ERROR: tried to set scale to a negative number");
                }
                HUD_SCALE = instruction->double_args[0];
            }
            break;
        case PRINT:
            {
                console_log(instruction->string_args[0]);
            }
            break;
        case FLIPPER: //layer, shaft x, shaft y, angle, angular range, length, large end radius, small end radius, ID out
            if (table)
            {
                //TODO (Phase VII) examine why backhands don't work well (may be realistic or may require higher friction on flipper)
                unsigned int layer = instruction->int_args[0];
                if (layer >= table->layers.size())
                {
                    console_log("SCRIPT ERROR: layer out of range");
                    return -1;
                }
                double shaft_x = instruction->double_args[0];
                double shaft_y = instruction->double_args[1];
                double angle = instruction->double_args[2];
                double angular_range = instruction->double_args[3];
                double length = instruction->double_args[4];
                double large_radius = instruction->double_args[5];
                double small_radius = instruction->double_args[6];
                if (large_radius < small_radius)
                {
                    double temp = small_radius;
                    small_radius = large_radius;
                    large_radius = temp;
                }
                double D = length - large_radius - small_radius;
                double Y = std::sqrt(D * D - (large_radius - small_radius) * (large_radius - small_radius) + small_radius * small_radius);
                double theta = std::acos((large_radius * large_radius + D * D - Y * Y) / (2.0 * large_radius * D));
                Point small_center = Point(shaft_x + D * std::cos(angle), shaft_y + D * std::sin(angle));
                Point p1 = Point(shaft_x + large_radius * std::cos(theta + angle), shaft_y + large_radius * std::sin(theta + angle));
                Point p2 = Point(small_center.x + small_radius * std::cos(theta + angle), small_center.y + small_radius * std::sin(theta + angle));
                Point p3 = Point(shaft_x + large_radius * std::cos(-theta + angle), shaft_y + large_radius * std::sin(-theta + angle));
                Point p4 = Point(small_center.x + small_radius * std::cos(-theta + angle), small_center.y + small_radius * std::sin(-theta + angle));
                ArcCollider* large_end = new ArcCollider(shaft_x, shaft_y, large_radius, angle + theta, angle + 2.0 * PI - theta);
                ArcCollider* small_end = new ArcCollider(small_center.x, small_center.y, small_radius, angle - theta, angle + theta);
                LineCollider* line1 = new LineCollider(p1.x, p1.y, p2.x, p2.y);
                LineCollider* line2 = new LineCollider(p4.x, p4.y, p3.x, p3.y);
                colliders.push_back(large_end);
                colliders.back()->id = colliders.size() - 1;
                colliders.push_back(small_end);
                colliders.back()->id = colliders.size() - 1;
                colliders.push_back(line1);
                colliders.back()->id = colliders.size() - 1;
                colliders.push_back(line2);
                colliders.back()->id = colliders.size() - 1;
                CompositeCollider* flipper_composite = new CompositeCollider;
                flipper_composite->add_collider(large_end);
                flipper_composite->add_collider(small_end);
                flipper_composite->add_collider(line1);
                flipper_composite->add_collider(line2);
                flipper_composite->add_endpoints();
                colliders.push_back(flipper_composite);
                colliders.back()->id = colliders.size() - 1;
                MobileCollider* flipper = new MobileCollider(flipper_composite, shaft_x, shaft_y, 0.0, 0.0, 0.0, 0.0, angular_range > 0.0 ? angular_range : 0.0, angular_range < 0.0 ? -angular_range : 0.0);
                flipper->angular_motion_bounce_coefficient = 1.25;
                flipper->live_catch_enabled = true;
                flipper->bounce_coefficient = 1.47;
                flipper->friction_coefficient = 1.92 * GLOBAL_FRICTION - 1.0;
                flipper->solid = true;
                colliders.push_back(flipper);
                colliders.back()->id = colliders.size() - 1;
                table->layers[layer]->add_collider(flipper);
                //TODO (Phase VII) default impact sound
                //TODO (Phase VII) default rotation stop sound
                int index = instruction->int_args[1];
                if (index != -1)
                {
                    integer_variables[index] = colliders.size() - 1;
                }
            } else
            {
                console_log("SCRIPT ERROR: no table declared");
                return -1;
            }
            break;
        default:
            console_log("SCRIPT WARNING: undefined command");
            return -1;
            break;
        }
        return 0;
    }
    void execute_from(unsigned int line)
    {
        //TODO (Phase X) test WHILE
        bool done = false;
        std::vector<Command> flow_stack;
        std::vector<int> gosub_stack;
        while (!done && line < instructions.size())
        {
            Instruction* instruction = instructions[line];
            if (execute(instruction) == -1)
            {
                std::string error_message = "Command #";
                error_message.append(std::to_string(line));
                console_log(error_message);
            }
            switch (instruction->command)
            {
                case IF:
                    if (instruction->int_args[0])
                    {
                        flow_stack.push_back(IF);
                    } else
                    {
                        bool skipping = true;
                        int nest = 0;
                        while (skipping)
                        {
                            ++line;
                            if (line >= instructions.size())
                            {
                                console_log("SCRIPT FLOW ERROR: IF without ENDIF");
                                return;
                            } else if (instructions[line]->command == IF)
                            {
                                ++nest;
                            } else if (instructions[line]->command == ELSE && !nest)
                            {
                                flow_stack.push_back(ELSE);
                                skipping = false;
                            } else if (instructions[line]->command == ENDIF)
                            {
                                if (nest)
                                {
                                    --nest;
                                } else
                                {
                                    skipping = false;
                                }
                            }
                        }
                    }
                    break;
                case ELSE:
                    if (flow_stack.back() == IF)
                    {
                        bool skipping = true;
                        int nest = 0;
                        while (skipping)
                        {
                            ++line;
                            if (line >= instructions.size())
                            {
                                console_log("SCRIPT FLOW ERROR: ELSE without ENDIF");
                                return;
                            } else if (instructions[line]->command == IF)
                            {
                                ++nest;
                            } else if (instructions[line]->command == ENDIF)
                            {
                                if (nest)
                                {
                                    --nest;
                                } else
                                {
                                    flow_stack.pop_back();
                                    skipping = false;
                                }
                            }
                        }
                    } else
                    {
                        console_log("SCRIPT FLOW ERROR: ELSE without IF");
                        std::string error_message = "Command #";
                        error_message.append(std::to_string(line));
                        console_log(error_message);
                        return;
                    }
                    break;
                case ENDIF:
                    if (flow_stack.back() == IF || flow_stack.back() == ELSE)
                    {
                        flow_stack.pop_back();
                    } else
                    {
                        console_log("SCRIPT FLOW ERROR: ENDIF without IF");
                        std::string error_message = "Command #";
                        error_message.append(std::to_string(line));
                        console_log(error_message);
                        return;
                    }
                    break;
                case WHILE:
                    if (instruction->int_args[0])
                    {
                        flow_stack.push_back(WHILE);
                    } else
                    {
                        bool skipping = true;
                        int nest = 0;
                        while (skipping)
                        {
                            ++line;
                            if (line >= instructions.size())
                            {
                                console_log("SCRIPT FLOW ERROR: WHILE without WEND");
                                return;
                            } else if (instructions[line]->command == WHILE)
                            {
                                ++nest;
                            } else if (instructions[line]->command == WEND)
                            {
                                if (nest)
                                {
                                    --nest;
                                } else
                                {
                                    skipping = false;
                                }
                            }
                        }
                    }
                    break;
                case WEND:
                    if (flow_stack.back() == WHILE)
                    {
                        flow_stack.pop_back();
                        int skipping = true;
                        while (skipping)
                        {
                            if (line < 0)
                            {
                                console_log("SCRIPT FLOW ERROR: WEND without WHILE");
                                std::string error_message = "Command #";
                                error_message.append(std::to_string(line));
                                console_log(error_message);
                                return;
                            } else if (instructions[line]->command == WHILE)
                            {
                                skipping = false;
                            }
                            --line;
                        }
                    } else
                    {
                        console_log("SCRIPT FLOW ERROR: WEND without WHILE");
                        std::string error_message = "Command #";
                        error_message.append(std::to_string(line));
                        console_log(error_message);
                        return;
                    }
                    break;
                case BREAK:
                    //TODO (Phase X) test breaking while loops
                    while (flow_stack.size() && flow_stack.back() != WHILE)
                    {
                        flow_stack.pop_back();
                    }
                    if (flow_stack.size() == 0)
                    {
                        done = true;
                    } else if (flow_stack.back() == WHILE) //TODO remove this condition (satisfied by previous loop) and test
                    {
                        bool skipping = true;
                        int nest = 0;
                        while (skipping)
                        {
                            ++line;
                            if (line >= instructions.size())
                            {
                                console_log("SCRIPT FLOW ERROR: WHILE without WEND");
                                return;
                            } else if (instructions[line]->command == WHILE)
                            {
                                ++nest;
                            } else if (instructions[line]->command == WEND)
                            {
                                if (nest)
                                {
                                    --nest;
                                } else
                                {
                                    skipping = false;
                                }
                            }
                        }
                    }
                    break;
                case GOSUB:
                    //TODO (Phase X) test (including nesting)
                    gosub_stack.push_back(line);
                    line = instruction->int_args[0];
                    flow_stack.push_back(GOSUB);
                    break;
                case RETURN:
                    while (flow_stack.size() && flow_stack.back() != GOSUB)
                    {
                        flow_stack.pop_back();
                    }
                    if (flow_stack.size() && flow_stack.back() == GOSUB)
                    {
                        flow_stack.pop_back();
                        line = gosub_stack.back();
                        gosub_stack.pop_back();
                    } else
                    {
                        console_log("SCRIPT FLOW ERROR: RETURN without GOSUB");
                        std::string error_message = "Command #";
                        error_message.append(std::to_string(line));
                        console_log(error_message);
                        return;
                    }
                    break;
                default:
                    break;
            }
            ++line;
        }
    }
    void pump_triggers()
    {
        if (table)
        {
            while (table->trigger_queue.size())
            {
                trigger_queue.push_back(table->trigger_queue[0]);
                table->trigger_queue.erase(table->trigger_queue.begin());
            }
        }
    }
    void pump_inputs()
    {
        if (input_handler)
        {
            while (input_handler->trigger_queue.size())
            {
                trigger_queue.push_back(input_handler->trigger_queue[0]);
                input_handler->trigger_queue.erase(input_handler->trigger_queue.begin());
            }
        }
    }
    void resolve_queue()
    {
        while (trigger_queue.size())
        {
            Event event = trigger_queue[0].first;
            int parameter = trigger_queue[0].second;
            if (event == COLLISION && !colliders[parameter]->solid)
            {
                //TODO look into preventing double triggers for non-drop target collisions
                //TODO replace erase with pop_front where applicable
                trigger_queue.pop_front();
                break;
            }
            if (trigger_points.find(event) != trigger_points.end())
            {
                std::unordered_map<int, int>* secondary_map = trigger_points[event];
                if (secondary_map->find(parameter) != secondary_map->end())
                {
                    int line = (*secondary_map)[parameter];
                    execute_from(line);
                }
            }
            trigger_queue.pop_front();
        }
    }
    void handle_resize()
    {
        for (unsigned int i = 0; i < backglasses.size(); ++i)
        {
            for (unsigned int j = 0; j < backglasses[i]->dot_matrix_displays.size(); ++j)
            {
                backglasses[i]->dot_matrix_displays[j].first->generate_alpha_mask(renderer);
            }
        }
    }
    void advance_timers()
    {
        //TODO (Phase X) test exact timing against off-by-one errors
        for (auto i = timers.begin(); i != timers.end(); ++i)
        {
            if (i->second)
            {
                int index = i->first;
                --integer_variables[index];
                if (integer_variables[index] == 0)
                {
                    trigger_queue.push_back(std::pair<Event, int>(TIMER_DONE, index));
                    timers[index] = false;
                }
            }
        }
    }
    void simulate_camera(double ms)
    {
        //TODO (Phase VIII) integrate analog camera control
        //TODO (Phase VIII) add backglass view bind
        if (input_handler->input_state[INPUT_CAMERA_HOME].first)
        {
            view->transition_time = 0.0;
            if (table)
            {
                view->center.x = table->dimensions.x / 2.0;
                view->center.y = table->dimensions.y / 2.0;
            } else
            {
                view->center.x = 0.0;
                view->center.y = 0.0;
            }
        }
        double dx = input_handler->input_state[ANALOG_CAMERA_X].second * CAMERA_PAN_SENSITIVITY * ms / 1000.0 / view->zoom;
        double dy = input_handler->input_state[ANALOG_CAMERA_Y].second * CAMERA_PAN_SENSITIVITY * ms / 1000.0 / view->zoom;
        double dz = (input_handler->input_state[INPUT_CAMERA_ZOOM_IN].first - input_handler->input_state[INPUT_CAMERA_ZOOM_OUT].first) * CAMERA_ZOOM_SENSITIVITY * ms / 1000.0;
        if (dx || dy || dz)
        {
            view->transition_time = 0.0;
        }
        view->center.x += dx;
        view->center.y += dy;
        view->zoom = std::min(std::max(view->zoom * std::pow(2.0, dz), 0.125), 1.0);
    }
};

void save_config(const char filename[], Script* script)
{
    //TODO (Phase X) test
    std::ofstream file(filename);
    file << DISPLAY_WIDTH << std::endl;
    file << DISPLAY_HEIGHT << std::endl;
    file << FRAMES_PER_SECOND << std::endl;
    file << FULLSCREEN << std::endl;
    file << BORDERLESS << std::endl;
    file << DISABLE_FLASHING << std::endl;
    file << std::to_string(UI_SCALE) << std::endl;
    file << std::to_string(HUD_SCALE) << std::endl;
    file << script->view->debug << std::endl;
    file << script->view->sprites << std::endl;
    file << script->view->HUD << std::endl;
    file << script->volume << std::endl;
    file << script->music_volume << std::endl;
    file << script->sound_volume << std::endl;
    file << script->mute << std::endl;
    file << CONTROLS_PATH << std::endl;
}

void load_config(const char filename[], Script* script)
{
    std::ifstream file(filename);
    std::string line;
    Instruction instruction;
    std::getline(file, line);
    instruction.command = RESOLUTION;
    instruction.int_args.push_back(std::stoi(line));
    std::getline(file, line);
    instruction.int_args.push_back(std::stoi(line));
    script->execute(&instruction);

    std::getline(file, line);
    int framerate = std::stoi(line);
    //TODO (Phase X) test various framerates
    instruction.command = FRAMERATE;
    instruction.int_args.clear();
    instruction.int_args.push_back(framerate);
    script->execute(&instruction);

    instruction.command = DISPLAY_MODE;
    instruction.int_args.clear();
    std::getline(file, line);
    instruction.int_args.push_back(std::stoi(line));
    std::getline(file, line);
    instruction.int_args.push_back(std::stoi(line));
    script->execute(&instruction);

    std::getline(file, line);
    DISABLE_FLASHING = std::stoi(line);

    std::getline(file, line);
    UI_SCALE = std::stod(line);

    std::getline(file, line);
    HUD_SCALE = std::stod(line);

    std::getline(file, line);
    script->view->debug = std::stoi(line);

    std::getline(file, line);
    script->view->sprites = std::stoi(line);

    std::getline(file, line);
    script->view->HUD = std::stoi(line);

    std::getline(file, line);
    script->volume = std::stoi(line);

    std::getline(file, line);
    script->music_volume = std::stoi(line);

    std::getline(file, line);
    script->sound_volume = std::stoi(line);

    std::getline(file, line);
    script->mute = std::stoi(line);

    std::getline(file, line);
    instruction.command = LOAD_CONTROLS;
    instruction.int_args.clear();
    instruction.string_args.push_back(CONTROLS_PATH);
    script->execute(&instruction);
}

void default_config(Script* script)
{
    Instruction instruction;
    instruction.command = RESOLUTION;
    instruction.int_args.push_back(640);
    instruction.int_args.push_back(480);
    script->execute(&instruction);

    int framerate = 60;
    FRAMES_PER_SECOND = framerate;
    MS_PER_FRAME = 1000.0 / framerate;
    PHYSICS_TICKS_PER_FRAME = 24 * 60 / framerate;

    instruction.command = DISPLAY_MODE;
    instruction.int_args.clear();
    instruction.int_args.push_back(0);
    instruction.int_args.push_back(0);
    script->execute(&instruction);

    DISABLE_FLASHING = false;

    UI_SCALE = 1.0;

    HUD_SCALE = 1.0;

    script->view->debug = false;

    script->view->sprites = true;

    script->view->HUD = true;

    script->volume = 128;

    script->music_volume = 128;

    script->sound_volume = 128;

    script->mute = false;

    instruction.command = LOAD_CONTROLS;
    instruction.int_args.clear();
    instruction.string_args.push_back("controls/controls.txt");
    script->execute(&instruction);
}

int main(int argc, char* argv[])
{
    //TODO investigate cases where ball triggers one-way hole collider without becoming trapped in hole
    uint64_t initial_count = SDL_GetPerformanceCounter();

    SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO | SDL_INIT_EVENTS | SDL_INIT_JOYSTICK);
    IMG_Init(IMG_INIT_PNG);
    Mix_OpenAudio(FREQUENCY, MIX_DEFAULT_FORMAT, OUTPUT_CHANNELS, SAMPLE_SIZE);
    SDL_Window* window = SDL_CreateWindow("Original Pinball Compilation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, DISPLAY_WIDTH, DISPLAY_HEIGHT, SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE | (SDL_WINDOW_FULLSCREEN * FULLSCREEN));
    //TODO change window title when loading table
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED); //TODO (Phase X) consider vsync
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "2");
    Mix_AllocateChannels(TRACK_COUNT * CHANNELS_PER_TRACK);
    SDL_GameController* gamepad = nullptr;
    if (SDL_NumJoysticks())
    {
        gamepad = SDL_GameControllerOpen(0);
    }
    for (int i = 0; i < TRACK_COUNT; ++i)
    {
        Mix_GroupChannels(i * CHANNELS_PER_TRACK, (i + 1) * CHANNELS_PER_TRACK - 1, i);
    }
    srand(time(NULL));

    bool quit = false;
    int frame_count = 0;
    View view;
    view.zoom = 0.25;

    Syntax syntax;
    InputHandler input_handler = InputHandler();
    Script script(renderer, window, &syntax, "scripts/Magical-Robot.txt");
    script.window = window;
    script.view = &view;
    script.input_handler = &input_handler;
    if (script.table)
    {
        view.center.x = script.table->dimensions.x / 2.0;
        view.center.y = script.table->dimensions.y / 2.0;
    }

    //TODO (Phase X) failsafe/handler for missing console font files
    FontSize console_font_size(renderer, "fonts/common/24high12wide.png", "fonts/common/24high12wide_kerningdata.txt"); //TODO (Phase VII) specialized large monospaced console font?
    Font console_font;
    console_font.add_size(&console_font_size);

    load_config("config/config.txt", &script);

    //For benchmarking
    uint64_t performance_frequency = SDL_GetPerformanceFrequency();
    uint64_t input_time = 0;
    uint64_t logic_time = 0;
    uint64_t rendering_time = 0;
    uint64_t delay_time = 0;
    uint64_t slowdown_time = 0;
    uint64_t load_time = SDL_GetPerformanceCounter() - initial_count;
    uint64_t script_time = 0;
    uint64_t physics_time = 0;

    while (!quit)
    {
        uint64_t frame_start = SDL_GetTicks();

        uint64_t performance_frame_start = SDL_GetPerformanceCounter();
        uint64_t input_time_this_frame = 0;
        uint64_t logic_time_this_frame = 0;
        uint64_t rendering_time_this_frame = 0;
        uint64_t delay_time_this_frame = 0;
        uint64_t slowdown_time_this_frame = 0;
        uint64_t script_time_this_frame = 0;
        uint64_t physics_time_this_frame = 0;

        SDL_Event e;
        while (SDL_PollEvent(&e))
        {
            //TODO (Phase X) optimize with switch/case statements where applicable
            if (e.type == SDL_QUIT)
            {
                quit = true;
            } else if (e.type == SDL_WINDOWEVENT && e.window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
            {
                DISPLAY_WIDTH = e.window.data1;
                DISPLAY_HEIGHT = e.window.data2;
                script.handle_resize();
            } else if (e.type == SDL_MOUSEMOTION && e.motion.state & SDL_BUTTON_LMASK)
            {
                script.view->transition_time = 0.0;
                script.view->center.x -= e.motion.xrel / script.view->zoom / PIXELS_PER_INCH;
                script.view->center.y -= e.motion.yrel / script.view->zoom / PIXELS_PER_INCH;
            } else if (e.type == SDL_MOUSEWHEEL)
            {
                script.view->transition_time = 0.0;
                script.view->zoom = std::min(std::max(script.view->zoom * std::pow(2.0, e.wheel.preciseY), 0.125), 1.0);
            } else if (e.type == SDL_JOYDEVICEADDED && e.cdevice.which == 0)
            {
                gamepad = SDL_GameControllerOpen(0);
            } else if (e.type == SDL_JOYDEVICEREMOVED && e.cdevice.which == 0)
            {
                //TODO (Phase X) test adding two controllers and unplugging first controller
                gamepad = nullptr;
            } else if (e.type == SDL_KEYDOWN && e.key.keysym.sym == SDLK_BACKQUOTE && !e.key.repeat)
            {
                script.show_console = !script.show_console;
                if (script.show_console)
                {
                    SDL_StartTextInput();
                } else
                {
                    SDL_StopTextInput();
                }
            } else if (script.show_console)
            {
                script.process_event(&e);
            } else
            {
                input_handler.process_event(&e);
            }
        }
        input_time_this_frame = SDL_GetPerformanceCounter() - performance_frame_start;

        for (int i = 0; i < PHYSICS_TICKS_PER_FRAME; ++i) //TODO examine memory leak (manifesting as DMD visual errors) caused when number of loops is not equal to PHYSICS_TICKS_PER_FRAME; see if time scale option is possible
        {
            uint64_t loop_start = SDL_GetPerformanceCounter();
            script.add_trigger(EVERY_TICK, 0);
            script.pump_triggers();
            script.pump_inputs();
            if (i == PHYSICS_TICKS_PER_FRAME - 1)
            {
                script.add_trigger(EVERY_FRAME, 0);
            }
            script.resolve_queue();
            script_time_this_frame += SDL_GetPerformanceCounter() - loop_start;
            uint64_t physics_start = SDL_GetPerformanceCounter();
            if (script.table)
            {
                script.table->simulate();
            }
            uint64_t physics_end = SDL_GetPerformanceCounter();
            physics_time_this_frame += physics_end - physics_start;
            script.pass_balls();
            script.advance_timers();
            script_time_this_frame += SDL_GetPerformanceCounter() - physics_end;
        }
        logic_time_this_frame = SDL_GetPerformanceCounter() - input_time_this_frame - performance_frame_start;

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        if (script.view)
        {
            script.view->pass_time(MS_PER_FRAME);
            script.simulate_camera(MS_PER_FRAME);
        }
        if (script.background)
        {
            SDL_Rect bg_crop;
            bg_crop.x = 0.0;
            bg_crop.y = 0.0;
            bg_crop.w = script.bg_width;
            bg_crop.h = script.bg_height;
            int center_x = bg_crop.w / 2;
            int center_y = bg_crop.h / 2;
            double bg_aspect_ratio = bg_crop.w * 1.0 / bg_crop.h;
            double display_aspect_ratio = DISPLAY_WIDTH * 1.0 / DISPLAY_HEIGHT;
            if (bg_aspect_ratio < display_aspect_ratio)
            {
                bg_crop.h = bg_crop.w / display_aspect_ratio;
            } else if (bg_aspect_ratio > display_aspect_ratio)
            {
                bg_crop.w = bg_crop.h * display_aspect_ratio;
            }
            bg_crop.x = center_x - bg_crop.w / 2;
            bg_crop.y = center_y - bg_crop.h / 2;
            SDL_RenderCopy(renderer, script.background, &bg_crop, nullptr);
        }
        if (script.table)
        {
            script.table->render(renderer, script.view);
        }
        if (script.show_console)
        {
            script.render_console(renderer, &console_font, 0, 0, 32);
        }
        SDL_RenderPresent(renderer);
        rendering_time_this_frame = SDL_GetPerformanceCounter() - logic_time_this_frame - input_time_this_frame - performance_frame_start;

        int now = SDL_GetTicks();
        double frame_end = frame_start + MS_PER_FRAME;
        if (now < frame_end)
        {
            delay_time_this_frame = (frame_end - now) * performance_frequency / 1000.0;
            SDL_Delay(frame_end - now);
        } else
        {
            slowdown_time_this_frame = (now - frame_end) * performance_frequency / 1000.0;
            std::cout << "WARNING: slowdown" << std::endl;
        }
        ++frame_count;

        input_time += input_time_this_frame;
        logic_time += logic_time_this_frame;
        rendering_time += rendering_time_this_frame;
        delay_time += delay_time_this_frame;
        slowdown_time += slowdown_time_this_frame;
        script_time += script_time_this_frame;
        physics_time += physics_time_this_frame;
    }
    uint64_t runtime = SDL_GetPerformanceCounter() - load_time - initial_count;

    Mix_Quit();
    IMG_Quit();
    SDL_Quit();
    if (gamepad)
    {
        SDL_GameControllerClose(gamepad);
    }

    std::cout << "Total load time          : " << load_time * 1000.0 / performance_frequency << " ms" << std::endl;
    std::cout << "Total runtime            : " << runtime * 1000.0 / performance_frequency << " ms" << std::endl;
    std::cout << std::endl;
    std::cout << "Time handling inputs     : " << input_time * 1000.0 / performance_frequency << " ms (" << input_time * 1000.0 / performance_frequency / frame_count << " per frame, " << input_time * 100.0 / runtime << "%)" << std::endl;
    std::cout << "Time computing game logic: " << logic_time * 1000.0 / performance_frequency << " ms (" << logic_time * 1000.0 / performance_frequency / frame_count << " per frame, " << logic_time * 100.0 / runtime << "%)" << std::endl;
    std::cout << "     Time running script : " << script_time * 1000.0 / performance_frequency << " ms (" << script_time * 1000.0 / performance_frequency / frame_count << " per frame, " << script_time * 100.0 / logic_time << "% of logic time)" << std::endl;
    std::cout << "     Time running physics: " << physics_time * 1000.0 / performance_frequency << " ms (" << physics_time * 1000.0 / performance_frequency / frame_count << " per frame, " << physics_time * 100.0 / logic_time << "% of logic time)" << std::endl;
    std::cout << "Time rendering           : " << rendering_time * 1000.0 / performance_frequency << " ms (" << rendering_time * 1000.0 / performance_frequency / frame_count << " per frame, " << rendering_time * 100.0 / runtime << "%)" << std::endl;
    std::cout << "Time waiting             : " << delay_time * 1000.0 / performance_frequency << " ms (" << delay_time * 1000.0 / performance_frequency / frame_count << " per frame, " << delay_time * 100.0 / runtime << "%)" << std::endl;
    std::cout << std::endl;
    std::cout << "Total slowdown time      : " << slowdown_time * 1000.0 / performance_frequency << " ms (" << slowdown_time * 1000.0 / performance_frequency / frame_count << " per frame, " << slowdown_time * 100.0 / runtime << "%)" << std::endl;

    return 0;
}
