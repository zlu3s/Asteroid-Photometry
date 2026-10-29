create table if not exists images (
    id integer primary key autoincrement,
    file_path text not null,
    filter text not null,
    time text not null,
    asteroid_value float not null,
    calibrated_value float not null,
    ref_star_name text not null,

    foreign key (ref_star_name) references ref_stars(name)
    foreign key (filter) references filters(name)
);

create table if not exists ref_stars (
    id integer primary key autoincrement,
    name text not null
);

create table if not exists filters (
    id integer primary key autoincrement,
    name text not null,
    count float default 0
);
