//! Typed Pub/Sub Helpers
//!
//! Convenience functions for publishing and subscribing to typed messages
//! over [`AppMessage::Custom`]. Data is serialized with MessagePack.

use serde::{de::DeserializeOwned, Serialize};

use crate::message::{AppMessage, MessageBus};

/// Publish a typed message on a named topic.
///
/// The data is serialized to MessagePack and sent as `AppMessage::Custom`.
pub fn publish<T: Serialize>(bus: &mut MessageBus, topic: &str, data: &T) {
    let payload = rmp_serde::to_vec(data).expect("failed to serialize topic payload");
    bus.send(AppMessage::Custom {
        topic: topic.to_string(),
        payload,
    });
}

/// Try to extract a typed message from a `Custom` event matching `topic`.
///
/// Returns `Some(T)` if the message matches and deserializes successfully,
/// `None` otherwise.
pub fn subscribe<T: DeserializeOwned>(msg: &AppMessage, topic: &str) -> Option<T> {
    if let AppMessage::Custom {
        topic: msg_topic,
        payload,
    } = msg
    {
        if msg_topic == topic {
            return rmp_serde::from_slice(payload).ok();
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_publish_subscribe_roundtrip() {
        let mut bus = MessageBus::new();

        #[derive(serde::Serialize, serde::Deserialize, Debug, PartialEq)]
        struct MyEvent {
            value: i32,
        }

        publish(&mut bus, "my_topic", &MyEvent { value: 42 });

        let messages = bus.drain_outbox();
        assert_eq!(messages.len(), 1);

        let result: Option<MyEvent> = subscribe(&messages[0], "my_topic");
        assert_eq!(result, Some(MyEvent { value: 42 }));
    }

    #[test]
    fn test_subscribe_wrong_topic() {
        let msg = AppMessage::Custom {
            topic: "other".into(),
            payload: rmp_serde::to_vec(&42i32).unwrap(),
        };
        let result: Option<i32> = subscribe(&msg, "my_topic");
        assert_eq!(result, None);
    }

    #[test]
    fn test_subscribe_non_custom_message() {
        let msg = AppMessage::RequestRedraw;
        let result: Option<i32> = subscribe(&msg, "any");
        assert_eq!(result, None);
    }
}
