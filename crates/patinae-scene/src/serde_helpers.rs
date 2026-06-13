//! Serde helpers for types that don't natively support serde.

pub mod mat4_serde {
    use lin_alg::f32::Mat4;
    use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

    pub fn serialize<S>(mat: &Mat4, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        mat.data.serialize(serializer)
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<Mat4, D::Error>
    where
        D: Deserializer<'de>,
    {
        let data: [f32; 16] = Deserialize::deserialize(deserializer)?;
        Ok(Mat4 { data })
    }
}

pub mod opt_mat4_serde {
    use lin_alg::f32::Mat4;
    use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

    pub fn serialize<S>(mat: &Option<Mat4>, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match mat {
            Some(m) => m.data.serialize(serializer),
            None => serializer.serialize_none(),
        }
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<Option<Mat4>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let opt: Option<[f32; 16]> = Deserialize::deserialize(deserializer)?;
        Ok(opt.map(|data| Mat4 { data }))
    }
}

pub mod vec3_serde {
    use lin_alg::f32::Vec3;
    use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

    pub fn serialize<S>(v: &Vec3, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        [v.x, v.y, v.z].serialize(serializer)
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<Vec3, D::Error>
    where
        D: Deserializer<'de>,
    {
        let [x, y, z]: [f32; 3] = Deserialize::deserialize(deserializer)?;
        Ok(Vec3::new(x, y, z))
    }
}

pub mod opt_vec3_serde {
    use lin_alg::f32::Vec3;
    use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

    pub fn serialize<S>(v: &Option<Vec3>, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match v {
            Some(v) => [v.x, v.y, v.z].serialize(serializer),
            None => serializer.serialize_none(),
        }
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<Option<Vec3>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let opt: Option<[f32; 3]> = Deserialize::deserialize(deserializer)?;
        Ok(opt.map(|[x, y, z]| Vec3::new(x, y, z)))
    }
}
